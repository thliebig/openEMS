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

#ifndef PROCESSCURRENT_H
#define PROCESSCURRENT_H

#include "processintegral.h"

class ProcessCurrent : public ProcessIntegral
{
public:
	ProcessCurrent(Engine_Interface_Base* eng_if);
	virtual ~ProcessCurrent();

	virtual std::string GetProcessingName() const {return "current integration";}

	virtual std::string GetIntegralName(int row) const;

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	//! Integrate currents flowing through an area
	virtual double CalcIntegral();

protected:
};

#endif // PROCESSCURRENT_H
