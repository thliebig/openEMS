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

#include "operator_mpi.h"

#include "operator_sse_compressed.h"
#include "engine_sse_compressed.h"
#include "engine_mpi.h"
#include "tools/array_ops.h"

Operator_MPI* Operator_MPI::New()
{
	cout << "Create FDTD operator (compressed SSE + MPI)" << endl;
	Operator_MPI* op = new Operator_MPI();
	op->Init();
	return op;
}

Operator_MPI::Operator_MPI() : Operator_SSE_Compressed()
{
}

Operator_MPI::~Operator_MPI()
{
	Delete();
}

Engine* Operator_MPI::CreateEngine() const
{
	Engine_MPI* eng = Engine_MPI::New(this);
	return eng;
}

void Operator_MPI::Init()
{
	Operator_sse::Init();
}

void Operator_MPI::Delete()
{
}

void Operator_MPI::Reset()
{
	Delete();
	Operator_sse::Reset();
}
