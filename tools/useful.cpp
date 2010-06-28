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

#include "useful.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>

unsigned int CalcNyquistNum(double fmax, double dT)
{
	if (fmax==0) return UINT_MAX;
	if (dT==0) return 1;
	double T0 = 1/fmax;
	return floor(T0/2/dT);
}

