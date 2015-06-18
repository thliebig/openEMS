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

#ifndef PROCESSMODEMATCH_H
#define PROCESSMODEMATCH_H

#include "processintegral.h"

class CSFunctionParser;

//! Processing class to match a mode to a given analytic function and return the integral value.
/*!
  The analytric function has to be definied in drawing units.
  It will return the integral value and the mode purity as a secondary value.
  */
class ProcessModeMatch : public ProcessIntegral
{
public:
	ProcessModeMatch(Engine_Interface_Base* eng_if);
	virtual ~ProcessModeMatch();

	virtual std::string GetProcessingName() const;

	virtual std::string GetIntegralName(int row) const;

	virtual void InitProcess();
	virtual void Reset();

	//! Set the field type (0 electric field, 1 magnetic field)
	void SetFieldType(int type);
	//! Set the mode function in the given direction ny. For example: SetModeFunction(0,"cos(pi/1000*x)*sin(pi/500*y)");
	void SetModeFunction(int ny, std::string function);

	virtual int GetNumberOfIntegrals() const {return 2;}
	virtual double* CalcMultipleIntegrals();

protected:
	//normal direction of the mode plane
	int m_ny;

	int m_ModeFieldType;

	double GetField(int ny, const unsigned int pos[3]);
	double GetEField(int ny, const unsigned int pos[3]);
	double GetHField(int ny, const unsigned int pos[3]);

	std::string m_ModeFunction[3];
	CSFunctionParser* m_ModeParser[2];

	unsigned int m_numLines[2];
	double** m_ModeDist[2];
};

#endif // PROCESSMODEMATCH_H
