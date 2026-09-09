/*
* Copyright (C) 2026 openEMS contributors
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#include "operator_metal.h"
#include "engine_metal.h"

using std::cout;
using std::endl;

Operator_Metal* Operator_Metal::New()
{
	cout << "Create FDTD operator (Metal field updates)" << endl;
	Operator_Metal* op = new Operator_Metal();
	op->Init();
	return op;
}

Operator_Metal::Operator_Metal() : Operator_sse()
{
}

Engine* Operator_Metal::CreateEngine()
{
	m_Engine = Engine_Metal::New(this);
	return m_Engine;
}
