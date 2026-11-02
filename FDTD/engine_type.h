/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_TYPE_H
#define ENGINE_TYPE_H

//! The engine types, of the engine (Engine::GetType()) and of the engine an operator creates (Operator::SetEngineType())
/*!
  A base of Engine, so the types read Engine::BASIC, Engine::SSE and Engine::GPU;
  on its own so that operator.h, which engine.h includes, can use them too.
  */
struct Engine_Types
{
	enum EngineType
	{
		BASIC, SSE, GPU, UNKNOWN
	};

	//! The name of \a type in the messages, e.g. "Create FDTD operator (<name>)"
	static const char* GetTypeName(EngineType type)
	{
		switch (type)
		{
		case BASIC:
			return "basic + multi-threading";
		case SSE:
			return "compressed SSE + multi-threading";
		case GPU:
			return "GPU";
		default:
			return "unknown";
		}
	}
};

#endif // ENGINE_TYPE_H
