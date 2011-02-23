/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY{} without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef ENGINE_MPI_H
#define ENGINE_MPI_H

#include "operator_mpi.h"
#include "engine_sse_compressed.h"
#include "mpi.h"

class Engine_MPI : public Engine_SSE_Compressed
{
public:
    Engine_MPI();

	static Engine_MPI* New(const Operator_MPI* op);
	virtual ~Engine_MPI();

	virtual void Init();
	virtual void Reset();

	virtual bool IterateTS(unsigned int iterTS);

protected:
	Engine_MPI(const Operator_MPI* op);
	const Operator_MPI* m_Op_MPI;

	MPI_Status stat;
	MPI_Request Send_Request[3];
	MPI_Request Recv_Request[3];

	//field buffer for MPI transfer...
	unsigned int m_BufferSize[3];
	float* m_BufferUp[3];
	float* m_BufferDown[3];

	//! Transfer all tangential voltages at the upper bounds to the lower bounds of the neighbouring MPI-processes
	virtual void SendReceiveVoltages();
	//! Transfer all tangential currents at the lower bounds to the upper bounds of the neighbouring MPI-processes
	virtual void SendReceiveCurrents();
};

#endif // ENGINE_MPI_H
