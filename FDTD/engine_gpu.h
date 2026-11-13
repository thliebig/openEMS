/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#ifndef ENGINE_GPU_H
#define ENGINE_GPU_H

#include "engine.h"
#include "gpu_backend.h"
#include "engine_plan.h"

#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

class GPU_Backend;
class GPU_Extension;

//! FDTD engine running the main updates on a GPU backend
/*!
  The fields live on the device (see GPU_Backend). The base class arrays
  (volt_ptr, curr_ptr) serve as host mirror in the basic engine layout, read by
  the field processing through GetVolt/GetCurr. If the device shares its
  memory with the host, the host mirror is a view of the device memory and no
  copies are needed.

  Every extension runs on the device, through its device implementation,
  written once for every backend (GPU_CreateExtension(), gpu_ext.h); an
  extension without one is an error. The fields stay on the device for a whole
  IterateTS() call. Without shared memory, the host mirror is not updated at the
  end of a batch: GetVolt/GetCurr read the z-line of a value from the device
  (e.g. for a probe), or the whole field if a batch reads many lines (e.g. for a
  dump).

  A cylindrical multi-grid level gets the engine of its sub-grid with
  Engine::SetSubGrid(): a GPU engine on a backend sharing the work stream of
  this one (GPU_Backend::NewSubGridBackend()), coupled on the device by
  SubGrid_Coupling_GPU (gpu_ext.h).

  The engine type is Engine::GPU, which the extension dispatcher and
  Engine_Interface_FDTD handle like the basic engine, because that is the layout
  of the host mirror. The engine owns the backend and the device extensions;
  it deletes the sub-grid first, then the extensions, then the backend.
  */
class Engine_GPU : public Engine
{
public:
	//! Create the engine for \a op on \a backend (see GPU_Backend::New), which the engine takes over, or for a sub-grid on a backend sharing the work stream of \a parent
	/*!
	  \a modelBuilt is called once the engine no longer needs the coefficients of
	  \a op: when the cell model is built, before the fields are allocated (the
	  lower memory peak, as Engine_CPU::New()).
	  */
	static Engine_GPU* New(const Operator* op, GPU_Backend* backend, GPU_Backend* parent=NULL, const std::function<void()>& modelBuilt=nullptr);
	virtual ~Engine_GPU();

	//! Create the backend and the device extensions
	virtual void Init();
	//! Free the device extensions before the backend they belong to
	virtual void Reset();

	//! A batch of timesteps, the host mirror is brought up to date at its end
	virtual bool IterateTS(unsigned int iterTS);

	//! Temporal blocking (see IterateTS()): the plans of the candidate tiles, and which runs
	/*!
	  If the extensions run on the device and there is no sub-grid: parallelograms along x
	  (Engine_BlockedPlan::Parallelograms(): the tiles in order, no wedges, half
	  the launches of the trapezoids), all of the device on one tile, the stream
	  order in place of the barriers. The tile has to fit the last-level cache:
	  of its size (GPU_Backend::CacheBytes(), or OPENEMS_GPU_CACHE in MiB) tiles
	  of 1/2, 3/4 and all of it, else of 8 to 256 MiB; they and the flat sweep
	  are tried for OPENEMS_GPU_TRIAL seconds each (default 0.3) at the start,
	  the fastest is kept. With the fused step (see
	  PrepareFused(): the choice for a small cache) the flat sweep is two
	  candidates, the fused step and the two updates, also without temporal
	  blocking; a tile or the two updates free its memory. OPENEMS_GPU_TBLOCK=0
	  switches blocking off, =k:W[,k:W...] runs those tiles, tried against each
	  other and the flat sweep only with OPENEMS_GPU_TRIAL (Engine_TilesFromEnv(),
	  as OPENEMS_CPU_TBLOCK).
	  */
	void PrepareBlocking();
	//! Start the trials of PrepareBlocking() with the first candidate
	void StartTrials();
	//! The steps of the next block, up to \a left: until an extension needs a flat step (as Engine_ThreadPool::BlockLength())
	unsigned int BlockLength(unsigned int left, unsigned int kmax) const;
	//! Run the k timesteps of \a plan
	void RunBlock(const Engine_BlockedPlan& plan);
	//! One step of a blocked plan: the order of the flat half step, restricted to the box (as Engine_ThreadPool::BlockedStep())
	void BlockedStep(const Engine_Step& step, unsigned int ts);
	//! Count \a steps timesteps of the running trial, move on when it took long enough
	void TrialProgress(unsigned int steps);
	//! Run (\a run) or pause the clock of the running trial; it runs only inside IterateTS(), not over the host's work between the batches (probes, dumps), and pausing waits for the device
	void TrialClock(bool run);

	//! The fused step of the backend (GPU_Backend::PrepareFused()) for the flat timesteps, if every extension tells the nodes its hooks between the main updates touch
	/*!
	  Not with a sub-grid or the extensions on the host. The trials of
	  PrepareBlocking() decide whether it runs. OPENEMS_GPU_FUSED=0 keeps the
	  two main updates separate, =1 runs the fused step without the trial
	  against them. Called by PrepareBlocking(): \a skip, if not empty, is why
	  the fused step is not worth trying (the cache holds tiles of the full
	  depth); =1 runs it anyway.
	  */
	void PrepareFused(const std::string& skip);
	//! One flat timestep with the fused step: the hooks in the order of VoltageHalfStep() and CurrentHalfStep()
	void FusedStep();
	//! One voltage half-step: extensions and main update, then the sub-grid and the coupling of the two (see Engine::SetSubGrid())
	void VoltageHalfStep();
	//! One current half-step: extensions and main update, then the sub-grid and the coupling
	void CurrentHalfStep();
	//! Count the completed timestep, of the sub-grid too
	void NextTimestep();
	//! After a batch of timesteps: fill this grid from the sub-grid, update the host mirror for the field processing
	void FinishBatch();

	//! The device backend holding the fields
	GPU_Backend* GetBackend() const {return m_Backend;}

	//! Sums of the squared voltages and currents on the device (see GPU_Backend::CalcFastEnergy()); false if it has none
	virtual bool CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) const;

	//! Host mirror of the fields (basic engine layout)
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& HostVoltages() const {return *volt_ptr;}
	ArrayLib::ArrayNIJK<FDTD_FLOAT>& HostCurrents() const {return *curr_ptr;}

	//! Make the device results visible in the host mirror
	void VoltagesToHost();
	void CurrentsToHost();
	//! Update the host mirror if it is out of date (see GetVolt())
	void UpdateHostMirror();

	//! Copy of the current fields in the snapshot, see GPU_Backend::SnapshotFields(); false if not available
	virtual bool SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr);
	//! Wait until the snapshot can be read, may be called from another thread
	void WaitSnapshot() const;
	//! Whether the snapshots copy only the raw values the dumps read (see AddSnapshotGather()), instead of the fields
	virtual bool CanSnapshotGather() const;
	//! Add the raw values \a region of a dump to the snapshots, at \a offset of the voltages or currents (\a h_field)
	virtual bool AddSnapshotGather(bool h_field, const Engine_GatherRegion& region, size_t& offset);

	//! Sums of the weighted raw values \a region on the device, see GPU_Backend::AddFieldDFT(); -1 if not available
	int AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count);
	void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights) {m_Backend->AccumulateFieldDFT(id, weights);}
	bool ReadFieldDFT(int id, unsigned int k, std::vector<std::complex<float>>& sums) {return m_Backend->ReadFieldDFT(id, k, sums);}

	//! Field values, read from the device if the host mirror is out of date (see class description)
	virtual FDTD_FLOAT GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const;
	virtual FDTD_FLOAT GetVolt(unsigned int n, const unsigned int pos[3]) const {return GetVolt(n, pos[0], pos[1], pos[2]);}
	virtual FDTD_FLOAT GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const;
	virtual FDTD_FLOAT GetCurr(unsigned int n, const unsigned int pos[3]) const {return GetCurr(n, pos[0], pos[1], pos[2]);}

protected:
	Engine_GPU(const Operator* op, GPU_Backend* backend, GPU_Backend* parent);

	GPU_Backend* m_ParentBackend;
	GPU_Backend* m_Backend;
	bool m_SharedMemory; //!< the host mirror is the device memory
	bool m_SnapshotsUsed; //!< SnapshotFields() succeeded once
	std::vector<Engine_GatherRegion> m_SnapshotGather[2]; //!< raw values of the snapshots: from the voltages, the currents (see AddSnapshotGather())
	size_t m_SnapshotGatherSize[2];                     //!< their number of raw values
	bool m_SnapshotGatherSet;                           //!< m_SnapshotGather was passed to the backend
	bool m_FieldDFTUsed;                                //!< AddFieldDFT() succeeded once
	//! How far the field dumps may use the device (OPENEMS_GPU_DUMPS, for tests and bug reports):
	//! 0 (device): device sums and snapshots, 1 (snapshot): no device sums, 2 (host): neither
	int m_DumpPath;

	//! device implementations of the engine extensions, same order as m_Eng_exts
	std::vector<GPU_Extension*> m_GPU_exts;

	//! Out-of-date host mirror of one field (see GetVolt()): the z-lines read from the device in this batch
	struct StaleField
	{
		bool stale;
		bool full_this_batch;   //!< this batch read the whole field
		bool full_last_batch;   //!< the last batch read the whole field, so will this one
		std::unordered_map<size_t, std::vector<FDTD_FLOAT>> lines;
	};
	mutable StaleField m_StaleVolt, m_StaleCurr;
	//! Read a value of a stale field from the device, keeping its z-line for the rest of the batch
	FDTD_FLOAT ReadStale(bool currents, unsigned int n, unsigned int x, unsigned int y, unsigned int z) const;
	//! The device updated the field, the host mirror is out of date from here on
	void MarkStale(StaleField& field);

	void ClearGPUExtensions();
	//! All lines of the grid: the main updates of a flat timestep
	Engine_Box Domain() const {return {{0, 0, 0}, {numLines[0], numLines[1], numLines[2]}};}
	//! The hooks of the device extensions around \a update, the main update of the voltages or \a currents of \a box at timestep \a ts, in the order of Engine
	template <typename Update>
	void DeviceHalfStep(bool currents, const Engine_Box& box, unsigned int ts, Update update);

	std::vector< std::vector<Engine_BlockedPlan> > m_BlockCandidates;   //!< per candidate tile its plans, every k from the largest down to 2
	int m_Block = -1;                                                  //!< the running candidate, -1 the flat sweep
	bool m_Fused = false;                                              //!< the flat timesteps run FusedStep()
	bool m_FusedForced = false;                                        //!< OPENEMS_GPU_FUSED=1: no trial against the two updates
	bool m_FlatTrial = false;                                          //!< the trials try the fused step and the two updates as the flat sweep
	bool m_FlatSplit = false;                                          //!< the trial of the two updates runs
	bool m_BlockingPrepared = false;                                   //!< PrepareBlocking() ran (at the first IterateTS())
	std::function<void()> m_ModelBuilt;                                //!< see New(), called once
	double m_NodeBytes = 0;                                            //!< GPU_CellModel::NodeBytes(), for the tiles
	//! the trials: candidate m_Trial (the tiles, then the flat sweep, then the two updates with m_FlatTrial), its timesteps and seconds; -1 done
	int m_Trial = -1;
	double m_TrialSeconds = 0.3;
	unsigned int m_TrialSteps = 0;
	double m_TrialTime = 0;        //!< the seconds of the running trial so far, inside IterateTS() only
	double m_TrialMark = -1;       //!< when the clock was started, -1 paused
	bool m_TrialWarm = false;      //!< the first steps of the running trial are still to come (not timed)
	std::vector<double> m_TrialSpeed;
};

#endif // ENGINE_GPU_H
