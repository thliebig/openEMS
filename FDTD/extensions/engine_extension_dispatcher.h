/*
*	Copyright (C) 2024 Yifeng Li <tomli@tomli.me>
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_EXTENSION_DISPATCHER_H
#define ENGINE_EXTENSION_DISPATCHER_H

#include "FDTD/engine.h"
#include "FDTD/engine_sse.h"

// In openEMS, all extensions are subclasses from the abstract Engine_Extension
// to implement features like Engine_Extension::Apply2Voltages(). When an
// Engine_Extension is created, it accepts a pointer to Engine. Simulation data
// is taken from the Engine using virtual functions like Engine::GetVolt(). This
// allows each Engine_Extension to work on all engines without consider their
// underlying data format.
//
// But this architecture created a serious performance problem: access to every
// cell requires a virtual function call, and performance was extremely low.
// As a workaround, each extension has the same simulation code repeated many
// times using non-virtual subclass engine functions in a big if/else check,
// e.g. 2 different Apply2Voltages(), one calls Engine::GetVolt(), another calls
// Engine_sse::GetVolt(). This led to code duplication and poor readability.
//
// This header solves the problem by the macro ENG_DISPATCH. Now within each
// extension, features are first written as macros like:
//
//     Apply2VoltagesImpl<Engine>(Engine*)
//
// which are used internally by:
//
//     Apply2Voltages()
//
// by calling:
//
//     ENG_DISPATCH(Apply2VoltagesImpl)
//     ENG_DISPATCH_ARGS(Apply2VoltagesImpl, arg1, arg2, ...)
//
// This is not pretty but it significantly reduces code duplication.
//
// In the future, all Extensions should probably be eventually converted be
// templates.

#define ENG_DISPATCH(impl) \
	switch (m_Eng->GetType()) \
	{ \
	case Engine::SSE: \
		(this)->template impl<Engine_sse>((Engine_sse*) m_Eng); \
		break; \
	case Engine::BASIC: \
		(this)->template impl<Engine>((Engine*) m_Eng); \
		break; \
	default: \
		/* requires change here if a new engine is added. */ \
		throw std::runtime_error( \
		    "Engine_Extension_Dispatcher: unknown engine!" \
		); \
	}

#define ENG_DISPATCH_ARGS(impl, ...) \
	switch (m_Eng->GetType()) \
	{ \
	case Engine::SSE: \
		(this)->template impl<Engine_sse>((Engine_sse*) m_Eng, __VA_ARGS__); \
		break; \
	case Engine::BASIC: \
		(this)->template impl<Engine>((Engine*) m_Eng, __VA_ARGS__); \
		break; \
	default: \
		/* requires change here if a new engine is added. */ \
		throw std::runtime_error( \
		    "Engine_Extension_Dispatcher: unknown engine!" \
		); \
	}

#endif // ENGINE_EXTENSION_DISPATCHER_H
