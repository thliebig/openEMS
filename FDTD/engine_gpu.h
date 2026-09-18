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
#include "operator_gpu.h"

class GPU_Backend;

//! FDTD engine running the main updates on a GPU backend
/*!
  The fields live on the device (see GPU_Backend). The base class arrays
  (volt_ptr, curr_ptr) serve as host mirror in the basic engine layout, read by
  the field processing through GetVolt/GetCurr and by the engine extensions.

  Extensions whose operator extension is not IsGPUSave() run on the host mirror.
  The engine then uploads both fields before and downloads the updated field
  after each main update (host fallback, correct but slow). Note that the pre-
  and post-update hooks may write the fields too, despite the Engine_Extension
  documentation: the UPML in all four, the TF/SF and the cylinder extension in
  the post-update hooks. If all extensions are GPU-save, the fields
  stay on the device for a whole IterateTS() call and the host mirror is
  updated at its end.

  A GPU-save extension is called through the same hooks and must work on the
  device data (see GetBackend()) unless FieldsOnHost() is true, in which case it
  runs on the host mirror like any other extension.
  */
class Engine_GPU : public Engine
{
public:
	static Engine_GPU* New(const Operator_GPU* op);
	virtual ~Engine_GPU();

	virtual void Init();
	virtual void Reset();

	virtual bool IterateTS(unsigned int iterTS);

	//! The device backend holding the fields
	GPU_Backend* GetBackend() const {return m_Backend;}

	//! True if the extensions run on the host mirror, which is then kept in sync during IterateTS()
	bool FieldsOnHost() const {return m_FieldsOnHost;}

protected:
	Engine_GPU(const Operator_GPU* op);

	GPU_Backend* m_Backend;
	bool m_FieldsOnHost;
};

#endif // ENGINE_GPU_H
