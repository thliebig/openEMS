/*
*	Copyright (C) 2010 Sebastian Held (sebastian.held@gmx.de)
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENGINE_THREADPOOL_H
#define ENGINE_THREADPOOL_H

#include <vector>

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>
#if defined(_MSC_VER)
#include <intrin.h>   // _mm_pause(), __yield()
#endif

#include <boost/thread.hpp>

#include "engine_plan.h"

class Engine;

//! A barrier that spins briefly, then yields, then sleeps: the phases of the blocked plans and the
//! extension hooks wait at one many times per timestep, often for a few microseconds only, but the
//! hooks that run on one thread (and the work between two calls) keep the others waiting longer
class Engine_SpinBarrier
{
public:
	explicit Engine_SpinBarrier(unsigned int count) : m_Count(count), m_Waiting(0), m_Generation(0), m_Sleeping(0) {}

	void wait()
	{
		const unsigned int gen = m_Generation.load(std::memory_order_acquire);
		if (m_Waiting.fetch_add(1, std::memory_order_acq_rel)+1 == m_Count)
		{
			// the last one releases the others, and with them the writes of all
			m_Waiting.store(0, std::memory_order_relaxed);
			m_Generation.fetch_add(1, std::memory_order_seq_cst);
			if (m_Sleeping.load(std::memory_order_seq_cst)>0)
			{
				// a sleeper checks the generation under the lock: taking it orders the wakeup after its wait
				{ std::lock_guard<std::mutex> lock(m_Mutex); }
				m_Wake.notify_all();
			}
			return;
		}
		// spin, then yield, for a time, not a count: a pause takes ~10 to ~140 cycles depending on the CPU
		const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
		std::chrono::steady_clock::duration waited(0);
		for (unsigned int i=0; m_Generation.load(std::memory_order_acquire)==gen; ++i)
		{
			if ((i&31)==31)
				waited = std::chrono::steady_clock::now()-t0;
			if (waited < std::chrono::microseconds(SPIN_US))
				Pause();
			else if (waited < std::chrono::microseconds(YIELD_US))
				std::this_thread::yield();   // more threads than cores, or a longer wait
			else
			{
				// a long wait: sleep, no core burnt (seq_cst against the release above: one of the two sees the other)
				std::unique_lock<std::mutex> lock(m_Mutex);
				m_Sleeping.fetch_add(1, std::memory_order_seq_cst);
				while (m_Generation.load(std::memory_order_seq_cst)==gen)
					m_Wake.wait(lock);
				m_Sleeping.fetch_sub(1, std::memory_order_relaxed);
				return;
			}
		}
	}

protected:
	//! the spin, then the yields, in microseconds (the spin measured as 2000 pauses on a Ryzen 7 5800X, ~27 us;
	//! small models up to 2x). A sleep after 1 ms of yields cost small models up to 16 % (woken after every
	//! gap between two calls); after 10 ms a wakeup is ~1 % of the wait. Enums, not static members:
	//! std::chrono takes them by reference
	enum {SPIN_US = 30, YIELD_US = 10000};
	//! the CPU's hint for a spin-wait loop: less power, and the core's other thread runs on
	static void Pause()
	{
#if defined(__x86_64__) || defined(__i386__)
		__builtin_ia32_pause();
#elif defined(_M_X64) || defined(_M_IX86)
		_mm_pause();
#elif defined(__aarch64__) || defined(__arm__)
		__asm__ __volatile__("yield");
#elif defined(_M_ARM64) || defined(_M_ARM)
		__yield();
#endif
	}
	const unsigned int m_Count;
	std::atomic<unsigned int> m_Waiting;
	std::atomic<unsigned int> m_Generation;
	std::atomic<unsigned int> m_Sleeping;   //!< the threads asleep in m_Wake
	std::mutex m_Mutex;
	std::condition_variable m_Wake;
};

//! The worker threads of a CPU engine
/*!
  Each thread updates its boxes of a plan (Engine_Plan, the flat sweep split
  along Engine::GetSplitAxis()) with Engine::UpdateVoltages() and
  Engine::UpdateCurrents() and runs the thread-aware hooks of the engine
  extensions, with a barrier after the main updates and after each extension.

  With blocked plans of the engine (Engine::BuildBlockedPlans(), a single grid)
  whole blocks of timesteps run phase by phase, a barrier after each, the
  extensions through their box hooks: the plan with the largest k that fits
  the timesteps left, then the smaller ones; what is left runs flat.

  A half step (voltages, then currents) is done for the engine and then for its
  chain of sub-grids (Engine::GetSubGridEngine(), cylindrical multi-grid), each
  with its own plan, followed by the coupling of the grids, innermost
  first (Engine::SyncSubGridVoltages(), Engine::SyncSubGridCurrents()); after the
  timesteps the grids are filled from their sub-grids (Engine::InterpolateFromSubGrid()).
  */
class Engine_ThreadPool
{
	friend class Engine_ThreadPool_Worker;
public:
	//! Threads for \a eng; \a numThreads 0: Engine::DefaultThreads(), or, if that is 0, one and more while the speed increases, see NextInterval()
	Engine_ThreadPool(Engine* eng, unsigned int numThreads);
	~Engine_ThreadPool();

	//! Run \a iterTS timesteps and return when they are done
	void IterateTS(unsigned int iterTS);

	//! The speed of the last interval, to find the best number of threads (with 0 threads requested)
	void NextInterval(float curr_speed);

	unsigned int GetNumThreads() const {return m_numThreads;}

protected:
	void ChangeNumThreads(unsigned int numThreads);
	void StopThreads();

	//! the half steps of grid \a g by thread \a threadID
	void VoltageHalfStep(size_t g, unsigned int threadID);
	void CurrentHalfStep(size_t g, unsigned int threadID);
	//! the coupling of all grids after a half step, by the first thread
	void SyncSubGrids(unsigned int threadID, bool voltages);

	void DoPreVoltageUpdates(Engine* eng, int threadID);
	void DoPostVoltageUpdates(Engine* eng, int threadID);
	void Apply2Voltages(Engine* eng, int threadID);
	void DoPreCurrentUpdates(Engine* eng, int threadID);
	void DoPostCurrentUpdates(Engine* eng, int threadID);
	void Apply2Current(Engine* eng, int threadID);

	Engine* m_Eng;
	std::vector<Engine*> m_Grids;   //!< the engine and its sub-grids, outermost first
	//! per grid the plan of the threads; fewer threads in it than in the pool on a small grid
	std::vector<Engine_Plan> m_Plans;
	//! the candidates of blocked plans of the engine (Engine::BuildBlockedPlans()), each by decreasing k
	Engine_BlockedChoice m_Choice;
	//! the candidate the workers run, -1 the flat sweep; changed between two calls of IterateTS() only
	int m_Active;
	bool m_Blocked;   //!< m_Active>=0
	const std::vector<Engine_BlockedPlan>& BlockedPlans() const {return m_Choice.candidates[m_Active];}
	void SetActive(int active) {m_Active = active; m_Blocked = (active>=0);}
	//! the trials of the candidates and the flat sweep at the start (Engine_BlockedChoice::trialSeconds)
	bool m_Trying;
	int m_Trial;                       //!< -1 the warm-up (the start candidate), then the candidates, then the flat sweep
	int m_TrialDir;                    //!< the direction of the search along the candidates, +1 or -1, 0 at its end
	std::vector<double> m_TrialSpeed;  //!< timesteps per second of each candidate and the flat sweep
	double m_TrialSeconds;
	unsigned long long m_TrialTS;
	//! account \a iterTS timesteps in \a seconds to the trial running, on to the next when it is complete
	void Trial(unsigned int iterTS, double seconds);
	unsigned long long m_TS_blocked, m_TS_flat;   //!< the timesteps run blocked and flat, for the log
	//! one step of the blocked plan at timestep \a ts, with the box hooks of the extensions
	void BlockedStep(const Engine_Step& step, unsigned int ts);
	//! throws if \a eng cannot update \a box (Engine::CanUpdateBox())
	static void CheckPlanBox(const Engine* eng, const Engine_Box& box, bool currents);
	//! the timesteps of the next block from timestep \a ts on, at most \a left; below 2: a flat step
	unsigned int BlockLength(unsigned int ts, unsigned int left) const;
	boost::thread_group* m_thread_group;
	boost::barrier* m_startBarrier;
	boost::barrier* m_stopBarrier;
	Engine_SpinBarrier* m_IterateBarrier;
	volatile unsigned int m_iterTS;
	volatile bool m_stopThreads;
	unsigned int m_numThreads;       //!< number of worker threads
	unsigned int m_max_numThreads;   //!< max. number of worker threads
	bool m_opt_speed;
	float m_last_speed;
};

#endif // ENGINE_THREADPOOL_H
