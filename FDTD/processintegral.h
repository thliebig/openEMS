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

#ifndef PROCESSINTEGRAL_H
#define PROCESSINTEGRAL_H

#include "processing.h"

//! Abstract base class for integral parameter processing
class ProcessIntegral : public Processing
{
public:
	virtual ~ProcessIntegral();

	virtual void InitProcess();
	virtual void FlushData();

protected:
	ProcessIntegral(Operator* op, Engine* eng);

	vector<FDTD_FLOAT> TD_Values;

	vector<_Complex double> FD_Values;
};

#endif // PROCESSINTEGRAL_H
