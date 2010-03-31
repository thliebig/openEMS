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

#ifndef ARRAY_OPS_H
#define ARRAY_OPS_H

#include "../FDTD/operator.h"

FDTD_FLOAT*** Create3DArray(const unsigned int* numLines);
void Delete3DArray(FDTD_FLOAT*** array, const unsigned int* numLines);

FDTD_FLOAT**** Create_N_3DArray(const unsigned int* numLines);
void Delete_N_3DArray(FDTD_FLOAT**** array, const unsigned int* numLines);

void Dump_N_3DArray2File(ostream &file, FDTD_FLOAT**** array, const unsigned int* numLines);

#endif // ARRAY_OPS_H
