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

#ifndef GPU_BACKEND_REFERENCE_H
#define GPU_BACKEND_REFERENCE_H

#include "gpu_backend.h"

//! Reference GPU backend running on the CPU
/*!
  Keeps its own "device" copies of the fields and coefficients, separate from the
  host mirror of Engine_GPU, so that all synchronization paths are exercised. The
  updates follow Engine::UpdateVoltages/UpdateCurrents operation by operation, the
  results are bit-identical to the basic engine. Used until a real device backend
  exists and as the reference to validate one against.
  */
class GPU_Backend_Reference : public GPU_Backend
{
public:
	GPU_Backend_Reference();
	virtual ~GPU_Backend_Reference() {}

	virtual std::string GetName() const {return "reference (CPU)";}

	virtual bool Init(const Operator* op);

	virtual void UpdateVoltages();
	virtual void UpdateCurrents();

	virtual void DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt);
	virtual void DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr);
	virtual void UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt);
	virtual void UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr);

	virtual GPU_Backend* NewSubGridBackend() {return new GPU_Backend_Reference();}

protected:
	unsigned int numLines[3];

	// "device" memory
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_volt;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_curr;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_vv;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_vi;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_ii;
	ArrayLib::ArrayNIJK<FDTD_FLOAT> m_iv;
};

#endif // GPU_BACKEND_REFERENCE_H
