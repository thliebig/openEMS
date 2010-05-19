/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef ENGINE_SSE_COMPRESSED_H
#define ENGINE_SSE_COMPRESSED_H

#include "engine_sse.h"
#include "operator_sse_compressed.h"

class Engine_SSE_Compressed : public Engine_sse
{
public:
	static Engine_SSE_Compressed* New(const Operator_SSE_Compressed* op);
	virtual ~Engine_SSE_Compressed();

protected:
	Engine_SSE_Compressed(const Operator_SSE_Compressed* op);
	const Operator_SSE_Compressed* Op;

	virtual void UpdateVoltages(unsigned int startX, unsigned int numX);
	virtual void UpdateCurrents(unsigned int startX, unsigned int numX);
};

#endif // ENGINE_SSE_COMPRESSED_H
