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

#ifndef OPERATOR_MPI_H
#define OPERATOR_MPI_H

#include "operator_sse_compressed.h"

class Operator_MPI : public Operator_SSE_Compressed
{
public:
	//! Create a new operator
	static Operator_MPI* New();
	virtual ~Operator_MPI();

	virtual Engine* CreateEngine() const;

protected:
	Operator_MPI();
	virtual void Init();
	void Delete();
	virtual void Reset();
};

#endif // OPERATOR_MPI_H
