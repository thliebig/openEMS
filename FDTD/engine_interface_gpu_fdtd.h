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

#ifndef ENGINE_INTERFACE_GPU_FDTD_H
#define ENGINE_INTERFACE_GPU_FDTD_H

#include "engine_interface_fdtd.h"
#include "operator_gpu.h"

class Engine_GPU;

//! Engine interface of the GPU engine, see Engine_Interface_FDTD.
/*!
  The fields live on the device, so the processing hooks of Engine_Interface_Base
  mean something here: the host mirror has to be up to date before the fields are
  read (PrepareFieldAccess()), and a dump can be evaluated from a snapshot the
  device copies while it goes on time-stepping (TakeFieldSnapshot()).
  The mirror has the basic engine layout, which is what lets CreateFieldGather()
  evaluate the dumped nodes from the raw arrays.
  */
class Engine_Interface_GPU_FDTD : public Engine_Interface_FDTD
{
public:
	Engine_Interface_GPU_FDTD(Operator_GPU* op);
	virtual ~Engine_Interface_GPU_FDTD();

	virtual double CalcFastEnergy() const;
	virtual void PrepareFieldAccess();
	virtual bool TakeFieldSnapshot(unsigned int slot, const float* &volt, const float* &curr);
	virtual void WaitFieldSnapshot(unsigned int slot) const;
	virtual bool PrepareSnapshotGather(Engine_Field_Gather* gather);
	virtual Engine_Field_Gather* CreateFieldGather(bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]) const;

protected:
	//! NULL for engines with another storage, e.g. a cylindrical operator running on the CPU
	Engine_GPU* m_Eng_GPU;
};

#endif // ENGINE_INTERFACE_GPU_FDTD_H
