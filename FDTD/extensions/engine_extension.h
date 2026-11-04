/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_EXTENSION_H
#define ENGINE_EXTENSION_H

#define ENG_EXT_PRIO_DEFAULT 0 //default engine extension priority

// priority definitions for some important extensions
#define ENG_EXT_PRIO_STEADYSTATE		+2e6  //steady state extension priority
#define ENG_EXT_PRIO_UPML				+1e6  //unaxial pml extension priority
#define ENG_EXT_PRIO_CYLINDER			+1e5  //cylindrial extension priority
#define ENG_EXT_PRIO_TFSF				+5e4  //total-field/scattered-field extension priority
#define ENG_EXT_PRIO_EXCITATION			-1000 //excitation priority

#include <string>

#include "FDTD/engine_plan.h"

class Operator_Extension;
class Engine;

//! Abstract base-class for all engine extensions
class Engine_Extension
{
public:
	virtual ~Engine_Extension();

	virtual void SetNumberOfThreads(int nrThread);

	//! This method will be called __before__ the main engine does the usual voltage updates. This method may __not__ change the engine voltages!!!
	virtual void DoPreVoltageUpdates() {}
	virtual void DoPreVoltageUpdates(int threadID);
	//! This method will be called __after__ the main engine does the usual voltage updates. This method may __not__ change the engine voltages!!!
	virtual void DoPostVoltageUpdates() {}
	virtual void DoPostVoltageUpdates(int threadID);
	//! This method will be called __after__ all updates to the voltages and extensions and may add/set its results to the engine voltages, but may __not__ rely on the current value of the engine voltages!!!
	virtual void Apply2Voltages() {}
	virtual void Apply2Voltages(int threadID);

	//! This method will be called __before__ the main engine does the usual current updates. This method may __not__ change the engine current!!!
	virtual void DoPreCurrentUpdates() {}
	virtual void DoPreCurrentUpdates(int threadID);
	//! This method will be called __after__ the main engine does the usual current updates. This method may __not__ change the engine current!!!
	virtual void DoPostCurrentUpdates() {}
	virtual void DoPostCurrentUpdates(int threadID);
	//! This method will be called __after__ all updates to the current and extensions and may add/set its results to the engine current, but may __not__ rely on the current value of the engine current!!!
	virtual void Apply2Current() {}
	virtual void Apply2Current(int threadID);

	//! Whether the extension runs in the blocked plans (Engine_BlockedPlan, temporal blocking; so far the CPU engine)
	/*!
	  Then the engine calls the box hooks below, per step on its box, instead of
	  the hooks above: a box hook touches only the cells in \a box and takes the
	  timestep from \a ts, not from Engine::GetNumberOfTimesteps(), which advances
	  by a whole block. An extension that does not support this keeps the engine
	  on the flat sweep.
	  */
	virtual bool SupportsBoxes() const {return false;}
	//! The bytes of the extension's per-cell state a timestep touches, for the tile size of the blocked plans
	virtual size_t BoxBytes() const {return 0;}
	//! Called once, single-threaded, before the first box hook, e.g. to index the extension's cells by box
	virtual void PrepareBoxes(const unsigned int numLines[3]) {(void)numLines;}
	//! Timestep \a ts must run flat, with the hooks above on a coherent grid (e.g. a whole-grid energy): blocks end before it
	virtual bool NeedsFlatStep(unsigned int ts) const {(void)ts; return false;}
	virtual void DoPreVoltageUpdatesBox(const Engine_Box& box, unsigned int ts) {(void)box; (void)ts;}
	virtual void DoPostVoltageUpdatesBox(const Engine_Box& box, unsigned int ts) {(void)box; (void)ts;}
	virtual void Apply2VoltagesBox(const Engine_Box& box, unsigned int ts) {(void)box; (void)ts;}
	virtual void DoPreCurrentUpdatesBox(const Engine_Box& box, unsigned int ts) {(void)box; (void)ts;}
	virtual void DoPostCurrentUpdatesBox(const Engine_Box& box, unsigned int ts) {(void)box; (void)ts;}
	virtual void Apply2CurrentBox(const Engine_Box& box, unsigned int ts) {(void)box; (void)ts;}

	//! Set the Engine to this extension. This will usually done automatically by Engine::AddExtension
	virtual void SetEngine(Engine* eng) {m_Eng=eng;}

	//! Get the operator extension this engine extension was created by
	Operator_Extension* GetOperatorExtension() const {return m_Op_ext;}

	//! Get the priority for this extension
	virtual int GetPriority() const {return m_Priority;}

	//! Set the priority for this extension
	virtual void SetPriority(int val) {m_Priority=val;}

	virtual bool operator< (const Engine_Extension& other);

	virtual std::string GetExtensionName() const;

protected:
	Engine_Extension(Operator_Extension* op_ext);

	Operator_Extension* m_Op_ext;
	Engine* m_Eng;

	int m_Priority;

	int m_NrThreads;
};

#endif // ENGINE_EXTENSION_H
