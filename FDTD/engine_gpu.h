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

#ifndef ENGINE_GPU_H
#define ENGINE_GPU_H

#include "engine.h"

#include <string>
#include <vector>

class GPU_Backend;
class GPU_Extension;

//! FDTD engine running the main updates on a GPU backend
/*!
  The fields live on the device (see GPU_Backend). The base class arrays
  (volt_ptr, curr_ptr) serve as host mirror in the basic engine layout, read by
  the field processing through GetVolt/GetCurr and by the engine extensions.
  If the device shares its memory with the host, the host mirror is a view of
  the device memory and no copies are needed.

  An extension runs on the device if its operator extension IsGPUSave() and the
  backend provides a device implementation (GPU_Backend::CreateExtension()).
  If all extensions do, the fields stay on the device for a whole IterateTS()
  call and the host mirror is updated at its end.

  Otherwise all extensions run on the host mirror, and the engine uploads both fields before and downloads the updated field
  after each main update (host fallback, correct but slow). Note that the pre-
  and post-update hooks may write the fields too, despite the Engine_Extension
  documentation: the UPML in all four, the TF/SF and the cylinder extension in
  the post-update hooks.
  */
class Engine_GPU : public Engine
{
public:
	//! Create the engine for \a op on the GPU backend \a backend (see GPU_Backend::New)
	static Engine_GPU* New(const Operator* op, const std::string& backend);
	virtual ~Engine_GPU();

	virtual void Init();
	virtual void Reset();

	virtual bool IterateTS(unsigned int iterTS);

	//! The device backend holding the fields
	GPU_Backend* GetBackend() const {return m_Backend;}

	//! True if the extensions run on the host mirror, which is then kept in sync during IterateTS()
	bool FieldsOnHost() const {return m_FieldsOnHost;}

protected:
	Engine_GPU(const Operator* op, const std::string& backend);

	std::string m_BackendName;
	GPU_Backend* m_Backend;
	bool m_FieldsOnHost;
	bool m_SharedMemory; //!< the host mirror is the device memory

	//! Make the device results visible in the host mirror
	void VoltagesToHost();
	void CurrentsToHost();
	//! Make the host mirror changes visible to the device
	void VoltagesToDevice();
	void CurrentsToDevice();

	//! device implementations of the engine extensions, same order as m_Eng_exts (fast path only)
	std::vector<GPU_Extension*> m_GPU_exts;

	void ClearGPUExtensions();
};

#endif // ENGINE_GPU_H
