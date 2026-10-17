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

#include <array>
#include <map>

#include "engine_interface_fdtd.h"

class Field_Gather_FDTD;

class Engine_GPU;

//! Engine interface of the GPU engine, see Engine_Interface_FDTD.
/*!
  The fields live on the device, so the processing hooks of Engine_Interface_Base
  mean something here: the host mirror has to be up to date before several threads
  read the fields (PrepareFieldAccess()), a dump can be evaluated from a snapshot
  the device takes while it goes on time-stepping (TakeFieldSnapshot()), and the
  device can keep the sums of a frequency domain dump (CreateFieldDFT()).
  The mirror has the basic engine layout, which is what lets CreateFieldGather()
  evaluate the dumped nodes from the raw arrays.

  It serves every operator whose engine is an Engine_GPU, with a Cartesian or a
  cylindrical mesh: what depends on the mesh is asked of the operator.
  */
class Engine_Interface_GPU_FDTD : public Engine_Interface_FDTD
{
public:
	Engine_Interface_GPU_FDTD(Operator* op);
	virtual ~Engine_Interface_GPU_FDTD();

	virtual double CalcFastEnergy() const;

	virtual void PrepareFieldAccess();
	virtual bool TakeFieldSnapshot(const float* &volt, const float* &curr);
	virtual void WaitFieldSnapshot() const;
	virtual bool PrepareSnapshotGather(Engine_Field_Gather* gather);
	virtual const void* GetEngineID() const {return m_Eng;}

	virtual int CreateFieldDFT(const Engine_Field_Gather* gather, unsigned int count);
	virtual void AccumulateFieldDFT(int id, const std::vector<std::complex<float>>& weights);
	virtual bool ReadFieldDFT(int id, std::vector<ArrayLib::ArrayNIJK<std::complex<float>>*>& fields);

	virtual Engine_Field_Gather* CreateFieldGather(bool h_field, const unsigned int numLines[3], unsigned int* const posLines[3]) const;

protected:
	Engine_GPU* m_Eng_GPU;

	//! the gathers of the field DFTs of the engine, by id (see CreateFieldDFT()); owned by their dumps
	std::map<int, const Field_Gather_FDTD*> m_FieldDFT;
};

#endif // ENGINE_INTERFACE_GPU_FDTD_H
