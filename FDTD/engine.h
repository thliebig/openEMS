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

#ifndef ENGINE_H
#define ENGINE_H

#include <fstream>
#include "operator.h"

namespace NS_Engine_Multithread
{
class thread; // evil hack to access numTS from multithreading context
}

class Engine_Extension;

class Engine
{
public:
	enum EngineType
	{
		BASIC, SSE, UNKNOWN
	};

	static Engine* New(const Operator* op);
	virtual ~Engine();

	virtual void Init();
	virtual void Reset();

	//!Iterate a number of timesteps
	virtual bool IterateTS(unsigned int iterTS);

	virtual unsigned int GetNumberOfTimesteps() {return numTS;}

	//this access functions muss be overloaded by any new engine using a different storage model
	inline virtual FDTD_FLOAT GetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z )		const { return volt[n][x][y][z]; }
	inline virtual FDTD_FLOAT GetVolt( unsigned int n, const unsigned int pos[3] )							const { return volt[n][pos[0]][pos[1]][pos[2]]; }
	inline virtual FDTD_FLOAT GetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z )		const { return curr[n][x][y][z]; }
	inline virtual FDTD_FLOAT GetCurr( unsigned int n, const unsigned int pos[3] )							const { return curr[n][pos[0]][pos[1]][pos[2]]; }

	inline virtual void SetVolt( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)	{ volt[n][x][y][z]=value; }
	inline virtual void SetVolt( unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )						{ volt[n][pos[0]][pos[1]][pos[2]]=value; }
	inline virtual void SetCurr( unsigned int n, unsigned int x, unsigned int y, unsigned int z, FDTD_FLOAT value)	{ curr[n][x][y][z]=value; }
	inline virtual void SetCurr( unsigned int n, const unsigned int pos[3], FDTD_FLOAT value )						{ curr[n][pos[0]][pos[1]][pos[2]]=value; }

	//! Execute Pre-Voltage extension updates
	virtual void DoPreVoltageUpdates();
	//! Main FDTD engine voltage updates
	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	//! Execute Post-Voltage extension updates
	virtual void DoPostVoltageUpdates();
	//! Apply extension voltage changes
	virtual void Apply2Voltages();

	//! Execute Pre-Current extension updates
	virtual void DoPreCurrentUpdates();
	//! Main FDTD engine current updates
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);
	//! Execute Post-Current extension updates
	virtual void DoPostCurrentUpdates();
	//! Apply extension current changes
	virtual void Apply2Current();

	inline size_t GetExtensionCount() {return m_Eng_exts.size();}
	inline Engine_Extension* GetExtension(size_t nr) {return m_Eng_exts.at(nr);}
	virtual void SortExtensionByPriority();

	EngineType GetType() const {return m_type;}

protected:
	EngineType m_type;

	Engine(const Operator* op);
	const Operator* Op;

	unsigned int numLines[3];

	FDTD_FLOAT**** volt;
	FDTD_FLOAT**** curr;
	unsigned int numTS;

	virtual void InitExtensions();
	virtual void ClearExtensions();
	vector<Engine_Extension*> m_Eng_exts;

	friend class NS_Engine_Multithread::thread; // evil hack to access numTS from multithreading context
};

#endif // ENGINE_H
