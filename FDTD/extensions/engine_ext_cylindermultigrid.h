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

#ifndef ENGINE_EXT_CYLINDERMULTIGRID_H
#define ENGINE_EXT_CYLINDERMULTIGRID_H

#include "FDTD/engine_cylindermultigrid.h"
#include "engine_extension.h"
#include "FDTD/operator_cylindermultigrid.h"

class Operator_Ext_CylinderMultiGrid;

class Engine_Ext_CylinderMultiGrid : public Engine_Extension
{
public:
	Engine_Ext_CylinderMultiGrid(Operator_Extension* op_ext, bool isBase);
	virtual ~Engine_Ext_CylinderMultiGrid();

	void SetBarrier(boost::barrier* waitBase, boost::barrier* waitChild, boost::barrier* waitSync);

	virtual void DoPreVoltageUpdates();
	virtual void DoPostVoltageUpdates();
	virtual void Apply2Voltages();

	virtual void DoPreCurrentUpdates();
	virtual void DoPostCurrentUpdates();
	virtual void Apply2Current();

	virtual void SetEngine(Engine* eng);

protected:
	void SyncVoltages();
	void SyncCurrents();

	Engine_CylinderMultiGrid* m_Eng_MG;

	boost::barrier *m_WaitOnBase;
	boost::barrier *m_WaitOnChild;
	boost::barrier *m_WaitOnSync;

	bool m_IsBase;
};

#endif // ENGINE_EXT_CYLINDERMULTIGRID_H
