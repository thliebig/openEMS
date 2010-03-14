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

#ifndef FDTD_EXAMPLES_H
#define FDTD_EXAMPLES_H

#include "ContinuousStructure.h"
#include "tinyxml.h"

void BuildDipol(const char* filename);

void BuildPlaneWave(const char* filename);

void BuildMSL(const char* filename);

void BuildCoaxial_Cartesian(const char* filename);

void BuildHelix(const char* filename);

#endif // FDTD_EXAMPLES_H
