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

#ifndef GPU_BACKEND_H
#define GPU_BACKEND_H

#include <string>

#include "tools/constants.h"
#include "tools/global.h"
#include "tools/arraylib/array_nijk.h"

class Operator;
class Engine;
class Engine_Extension;

//! Device implementation of an engine extension, see GPU_Backend::CreateExtension()
/*!
  Engine_GPU calls it through the same hooks and in the same order as the engine
  extension it replaces, while the fields stay on the device. It works on the
  device data of the backend that created it.
  */
class GPU_Extension
{
public:
	virtual ~GPU_Extension() {}

	virtual void DoPreVoltageUpdates() {}
	virtual void DoPostVoltageUpdates() {}
	virtual void Apply2Voltages() {}

	virtual void DoPreCurrentUpdates() {}
	virtual void DoPostCurrentUpdates() {}
	virtual void Apply2Current() {}

	//! Called after each IterateTS() of the engine, with the host mirror up to date
	virtual void Synchronize() {}
};

//! Abstract interface to the device used by Engine_GPU
/*!
  A backend owns the device memory: the voltage and current fields and the update
  coefficients of the operator. It performs the main FDTD updates on the device.

  Engine_GPU keeps a host mirror of the fields in the basic Engine layout (ArrayNIJK),
  used by the field processing and by extensions without a device implementation,
  and synchronizes it with the Download/Upload methods.
  */
class GPU_Backend
{
public:
	//! Create a backend by name: "reference" or "auto" (the best device backend available). It's the responsibility of the caller to free it.
	static GPU_Backend* New(const std::string& name);

	virtual ~GPU_Backend() {}

	virtual std::string GetName() const = 0;

	//! Allocate the device fields (all zero) and upload the update coefficients of \a op
	virtual bool Init(const Operator* op) = 0;

	//! Main voltage update over the whole domain, equivalent to Engine::UpdateVoltages
	virtual void UpdateVoltages() = 0;
	//! Main current update over the whole domain, equivalent to Engine::UpdateCurrents
	virtual void UpdateCurrents() = 0;

	//! Copy the device voltages to the host mirror
	virtual void DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt) = 0;
	//! Copy the device currents to the host mirror
	virtual void DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr) = 0;
	//! Copy the host mirror voltages to the device
	virtual void UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt) = 0;
	//! Copy the host mirror currents to the device
	virtual void UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr) = 0;

	//! Create the device implementation of the engine extension \a eng_ext of engine \a eng, or NULL if this backend has none
	virtual GPU_Extension* CreateExtension(Engine_Extension* eng_ext, Engine* eng) {UNUSED(eng_ext); UNUSED(eng); return NULL;}
};

#endif // GPU_BACKEND_H
