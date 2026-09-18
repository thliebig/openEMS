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

#include "engine_gpu.h"
#include "operator_gpu.h"

using std::cout;
using std::endl;

Operator_GPU* Operator_GPU::New()
{
	cout << "Create FDTD operator (GPU)" << endl;
	Operator_GPU* op = new Operator_GPU();
	op->Init();
	return op;
}

Operator_GPU::Operator_GPU() : Operator()
{
}

Operator_GPU::~Operator_GPU()
{
}

Engine* Operator_GPU::CreateEngine()
{
	m_Engine = Engine_GPU::New(this);
	return m_Engine;
}
