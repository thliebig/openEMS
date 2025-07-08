/*
*	Copyright (C) 2010 Sebastian Held (Sebastian.Held@gmx.de)
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

#include "engine_sse.h"
#include "operator_sse.h"
#include "tools/array_ops.h"
//#include "processfields.h"

Operator_sse* Operator_sse::New()
{
	cout << "Create FDTD operator (SSE)" << endl;
	Operator_sse* op = new Operator_sse();
	op->Init();
	return op;
}

Operator_sse::Operator_sse() : Operator()
{
	f4_vv_ptr = NULL;
	f4_vi_ptr = NULL;
	f4_iv_ptr = NULL;
	f4_ii_ptr = NULL;
}

Operator_sse::~Operator_sse()
{
	Delete();
}

Engine* Operator_sse::CreateEngine()
{
	//! create a special sse-engine
	m_Engine = Engine_sse::New(this);
	return m_Engine;
}

void Operator_sse::Init()
{
	Operator::Init();
	f4_vv_ptr = NULL;
	f4_vi_ptr = NULL;
	f4_iv_ptr = NULL;
	f4_ii_ptr = NULL;
}

void Operator_sse::Delete()
{
	delete f4_vv_ptr;
	delete f4_vi_ptr;
	delete f4_iv_ptr;
	delete f4_ii_ptr;
	f4_vv_ptr = NULL;
	f4_vi_ptr = NULL;
	f4_iv_ptr = NULL;
	f4_ii_ptr = NULL;
}

void Operator_sse::Reset()
{
	Delete();
	Operator::Reset();
}


void Operator_sse::InitOperator()
{
	delete f4_vv_ptr;
	delete f4_vi_ptr;
	delete f4_iv_ptr;
	delete f4_ii_ptr;

	numVectors =  ceil((double)numLines[2]/4.0);

	f4_vv_ptr = new ArrayLib::ArrayNIJK<f4vector>(
		"f4_vv", {numLines[0], numLines[1], numVectors}
	);
	f4_vi_ptr = new ArrayLib::ArrayNIJK<f4vector>(
		"f4_vi", {numLines[0], numLines[1], numVectors}
	);
	f4_iv_ptr = new ArrayLib::ArrayNIJK<f4vector>(
		"f4_iv", {numLines[0], numLines[1], numVectors}
	);
	f4_ii_ptr = new ArrayLib::ArrayNIJK<f4vector>(
		"f4_ii", {numLines[0], numLines[1], numVectors}
	);
}
