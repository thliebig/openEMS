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

#ifndef PROCESSVOLTAGE_H
#define PROCESSVOLTAGE_H

#include "processintegral.h"

//! Process voltage along a line from start to stop coordinates. ATM integration along the axis e.g.: in x, then y then z direction (Future: diagonal integration)
class ProcessVoltage : public ProcessIntegral
{
public:
	ProcessVoltage(Engine_Interface_Base* eng_if);
	virtual ~ProcessVoltage();

	virtual std::string GetProcessingName() const {return "voltage integration";}

	virtual std::string GetIntegralName(int row) const;

	virtual void DefineStartStopCoord(double* dstart, double* dstop);

	virtual double CalcIntegral();

protected:
};

#endif // PROCESSVOLTAGE_H
