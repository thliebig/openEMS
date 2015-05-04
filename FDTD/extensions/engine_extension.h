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
#define ENG_EXT_PRIO_CYLINDERMULTIGRID	-3000 //cylindrial multi-grid extension priority

#include <string>

class Operator_Extension;
class Engine;

//! Abstract base-class for all engine extensions
class Engine_Extension
{
public:
	virtual ~Engine_Extension();

	virtual void SetNumberOfThreads(int nrThread);

	//! This methode will be called __before__ the main engine does the usual voltage updates. This methode may __not__ change the engine voltages!!!
	virtual void DoPreVoltageUpdates() {}
	virtual void DoPreVoltageUpdates(int threadID);
	//! This methode will be called __after__ the main engine does the usual voltage updates. This methode may __not__ change the engine voltages!!!
	virtual void DoPostVoltageUpdates() {}
	virtual void DoPostVoltageUpdates(int threadID);
	//! This methode will be called __after__ all updates to the voltages and extensions and may add/set its results to the engine voltages, but may __not__ rely on the current value of the engine voltages!!!
	virtual void Apply2Voltages() {}
	virtual void Apply2Voltages(int threadID);

	//! This methode will be called __before__ the main engine does the usual current updates. This methode may __not__ change the engine current!!!
	virtual void DoPreCurrentUpdates() {}
	virtual void DoPreCurrentUpdates(int threadID);
	//! This methode will be called __after__ the main engine does the usual current updates. This methode may __not__ change the engine current!!!
	virtual void DoPostCurrentUpdates() {}
	virtual void DoPostCurrentUpdates(int threadID);
	//! This methode will be called __after__ all updates to the current and extensions and may add/set its results to the engine current, but may __not__ rely on the current value of the engine current!!!
	virtual void Apply2Current() {}
	virtual void Apply2Current(int threadID);

	//! Set the Engine to this extention. This will usually done automatically by Engine::AddExtension
	virtual void SetEngine(Engine* eng) {m_Eng=eng;}

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
