/*
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

#ifndef ENGINE_EXTENSION_PHASE_DISPATCH_H
#define ENGINE_EXTENSION_PHASE_DISPATCH_H

#include <array>
#include <cstddef>
#include <typeinfo>
#include <vector>

// Per-phase extension schedules for Engine_Multithread.
//
// The multithreaded engine calls every extension in every one of the six
// extension phases and waits at a barrier after each call, even when the
// extension's hook is empty. A schedule lists, per phase, only the extensions
// that have work to do. A run of skipped extensions keeps one barrier, so the
// order of the remaining hooks and the barriers between them is unchanged.
namespace EngineExtensionPhaseDispatch
{

enum Phase
{
	PRE_VOLTAGE = 0,
	POST_VOLTAGE,
	APPLY_VOLTAGE,
	PRE_CURRENT,
	POST_CURRENT,
	APPLY_CURRENT,
	PHASE_COUNT
};

typedef unsigned char PhaseMask;
static const PhaseMask ALL_PHASES = (1u << PHASE_COUNT) - 1u;

inline PhaseMask Bit(Phase phase)
{
	return static_cast<PhaseMask>(1u << phase);
}

//! Phases in which an extension has a non-empty hook.
//! Only the exact dynamic types are recognized. A derived or unknown
//! extension keeps all phases, so inherited hooks are never skipped.
template <typename Base, typename Upml, typename Excitation, typename LumpedRLC, typename MurABC>
PhaseMask ActivePhases(const Base* extension)
{
	if (typeid(*extension) == typeid(Upml))
		return Bit(PRE_VOLTAGE) | Bit(POST_VOLTAGE) | Bit(PRE_CURRENT) | Bit(POST_CURRENT);
	if (typeid(*extension) == typeid(Excitation))
		return Bit(APPLY_VOLTAGE) | Bit(APPLY_CURRENT);
	if (typeid(*extension) == typeid(LumpedRLC))
		return Bit(PRE_VOLTAGE) | Bit(APPLY_VOLTAGE);
	if (typeid(*extension) == typeid(MurABC))
		return Bit(PRE_VOLTAGE) | Bit(POST_VOLTAGE) | Bit(APPLY_VOLTAGE);
	return ALL_PHASES;
}

template <typename Extension>
struct PhaseLists
{
	// A NULL entry is the single barrier kept for a run of skipped
	// extensions. Any other entry is a hook call followed by its barrier.
	std::array<std::vector<Extension*>, PHASE_COUNT> phases;

	void clear()
	{
		for (std::size_t n = 0; n < phases.size(); ++n)
			phases[n].clear();
	}
};

//! The pre-update phases run the extensions in reverse priority order.
inline bool IsReversePriorityPhase(Phase phase)
{
	return phase == PRE_VOLTAGE || phase == PRE_CURRENT;
}

//! Build the schedules from the priority-sorted extension list.
template <typename Extension, typename Classifier>
void Build(const std::vector<Extension*>& sortedExtensions, PhaseLists<Extension>& lists, Classifier classify)
{
	lists.clear();
	for (int phaseIndex = 0; phaseIndex < PHASE_COUNT; ++phaseIndex)
	{
		const Phase phase = static_cast<Phase>(phaseIndex);
		const PhaseMask phaseBit = Bit(phase);
		std::vector<Extension*>& schedule = lists.phases[phaseIndex];
		bool pendingEmptyRun = false;
		const auto appendInOrder = [&](Extension* extension) {
			if (classify(extension) & phaseBit)
			{
				if (pendingEmptyRun)
				{
					schedule.push_back(static_cast<Extension*>(0));
					pendingEmptyRun = false;
				}
				schedule.push_back(extension);
			}
			else
				pendingEmptyRun = true;
		};
		if (IsReversePriorityPhase(phase))
		{
			for (std::size_t n = sortedExtensions.size(); n > 0; --n)
				appendInOrder(sortedExtensions[n - 1]);
		}
		else
		{
			for (std::size_t n = 0; n < sortedExtensions.size(); ++n)
				appendInOrder(sortedExtensions[n]);
		}
		if (pendingEmptyRun)
			schedule.push_back(static_cast<Extension*>(0));
	}
}

template <typename Extension, typename CallHook, typename WaitAtBarrier>
void Dispatch(const std::vector<Extension*>& schedule, CallHook callHook, WaitAtBarrier waitAtBarrier)
{
	for (std::size_t n = 0; n < schedule.size(); ++n)
	{
		if (schedule[n])
			callHook(schedule[n]);
		waitAtBarrier();
	}
}

} // namespace EngineExtensionPhaseDispatch

#endif // ENGINE_EXTENSION_PHASE_DISPATCH_H
