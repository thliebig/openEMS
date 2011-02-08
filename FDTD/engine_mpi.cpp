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

#include "engine_mpi.h"

Engine_MPI* Engine_MPI::New(const Operator_MPI* op)
{
	cout << "Create FDTD engine (compressed SSE + MPI)" << endl;
	Engine_MPI* e = new Engine_MPI(op);
	e->Init();
	return e;
}

Engine_MPI::Engine_MPI(const Operator_MPI* op) : Engine_SSE_Compressed(op)
{
	m_Op_MPI = op;
}

Engine_MPI::~Engine_MPI()
{
	Reset();
}

void Engine_MPI::Init()
{
	Engine_SSE_Compressed::Init();

}

void Engine_MPI::Reset()
{
	Engine_SSE_Compressed::Reset();
}
