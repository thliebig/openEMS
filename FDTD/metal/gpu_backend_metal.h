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

#ifndef GPU_BACKEND_METAL_H
#define GPU_BACKEND_METAL_H

#include <vector>

#include "FDTD/gpu_backend.h"

//! GPU backend using Apple Metal
/*!
  The fields and coefficients live in shared-storage Metal buffers (unified
  memory), which Engine_GPU uses directly as its host mirror. All kernels of a
  batch of timesteps are encoded into one serial compute encoder, which is only
  committed when the host needs the fields.

  This header is plain C++, the Metal state is kept in Impl (see metal_internal.h).
  */
class GPU_Backend_Metal : public GPU_Backend
{
public:
	//! Create the backend on the default Metal device, or return NULL if there is none
	static GPU_Backend_Metal* New();
	virtual ~GPU_Backend_Metal();

	virtual std::string GetName() const;

	virtual bool Init(const Operator* op);

	virtual void UpdateVoltages();
	virtual void UpdateCurrents();

	virtual void DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt);
	virtual void DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr);
	virtual void UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt);
	virtual void UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr);

	virtual FDTD_FLOAT* GetSharedVoltages() const;
	virtual FDTD_FLOAT* GetSharedCurrents() const;
	virtual void Synchronize();

	virtual GPU_Backend* NewSubGridBackend();
	virtual GPU_MultiGridLink* CreateMultiGridLink(GPU_Backend* sub_grid, const GPU_MultiGridInterpolation& interpol);

	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, Engine* eng);

	struct Impl;
	Impl* GetImpl() const {return d;}

protected:
	GPU_Backend_Metal(Impl* impl);

	//! Store the coefficients as distinct sets and a set index per node, returns false if there are too many sets
	bool CompressCoefficients(const std::vector<float>& vv, const std::vector<float>& vi,
	                          const std::vector<float>& ii, const std::vector<float>& iv);

	Impl* d;
};

#endif // GPU_BACKEND_METAL_H
