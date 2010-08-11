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

	//! Flush FD data to disk
	virtual void FlushData();

	//! This method should calculate the integral parameter and must be overloaded for each derived class
	virtual double CalcIntegral() {return 0;}

	//! This method will write the TD and FD dump files using CalcIntegral() to calculate the integral parameter
	virtual int Process();

protected:
	ProcessIntegral(Operator* op, Engine* eng);

	//! timeshift to be used in TD and FD data, e.g. 0.5*dT in case of current based parameter
	double m_TimeShift;

	vector<FDTD_FLOAT> TD_Values;
	vector<_Complex double> FD_Values;
};

#endif // PROCESSINTEGRAL_H
