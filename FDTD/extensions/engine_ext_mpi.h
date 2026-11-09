/*
*	Copyright (C) 2026 Ivan Iudice
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

#ifndef ENGINE_EXT_MPI_H
#define ENGINE_EXT_MPI_H

#include <memory>
#include <vector>

#include "engine_extension.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"

//! The lowest priority: the exchange sees the fields after all other extensions changed them
#define ENG_EXT_PRIO_MPI -1e7

class Operator_Ext_MPI;
class HaloThreads;

//! Copies the shared planes to the neighbour ranks after every half step (see MPI_Domain)
/*!
  The directions are exchanged one after the other, so the corner lines of a
  block reach the diagonal neighbours through two exchanges.
  */
class Engine_Ext_MPI : public Engine_Extension
{
public:
	Engine_Ext_MPI(Operator_Ext_MPI* op_ext);
	virtual ~Engine_Ext_MPI();

	virtual void Apply2Voltages();
	virtual void Apply2Current();

protected:
	Operator_Ext_MPI* m_Op_Ext_MPI;

	std::vector<float> m_SendBuf[3];
	std::vector<float> m_RecvBuf[3];

	//! Threads that copy the planes with the calling engine thread: the other engine threads wait meanwhile
	std::unique_ptr<HaloThreads> m_Threads;

	template <typename EngType>
	void ExchangeImpl(EngType* eng, bool voltages);

	//! Copy the two tangential components of the rows [\a r0, \a r1) of the plane \a pos (direction \a n) to/from \a buf
	template <typename EngType>
	void PackPlane(EngType* eng, bool voltages, int n, unsigned int pos, float* buf, unsigned int r0, unsigned int r1) const;
	template <typename EngType>
	void UnpackPlane(EngType* eng, bool voltages, int n, unsigned int pos, const float* buf, unsigned int r0, unsigned int r1) const;
	//! Copy the plane on the threads (pack or unpack)
	template <typename EngType>
	void CopyPlane(EngType* eng, bool voltages, bool pack, int n, unsigned int pos, float* buf);
};

#endif // ENGINE_EXT_MPI_H
