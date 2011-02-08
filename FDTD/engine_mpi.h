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

class Engine_MPI : public Engine_SSE_Compressed
{
public:
    Engine_MPI();

	static Engine_MPI* New(const Operator_MPI* op);
	virtual ~Engine_MPI();

	virtual void Init();
	virtual void Reset();

protected:
	Engine_MPI(const Operator_MPI* op);
	const Operator_MPI* m_Op_MPI;

};

#endif // ENGINE_MPI_H
