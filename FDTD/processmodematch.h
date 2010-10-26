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

class ProcessModeMatch : public ProcessIntegral
{
public:
	ProcessModeMatch(Operator* op, Engine* eng);
	virtual ~ProcessModeMatch();

	virtual void InitProcess();
	virtual void Reset();

	void SetFieldType(int type);
	void SetModeFunction(int ny, string function);

	virtual int GetNumberOfIntegrals() const {return 2;}
	virtual double* CalcMultipleIntegrals();

protected:
	//normal direction of the mode plane
	int m_ny;

	int m_ModeFieldType;

	double GetField(int ny, const unsigned int pos[3]);
	double GetEField(int ny, const unsigned int pos[3]);
	double GetHField(int ny, const unsigned int pos[3]);

	string m_ModeFunction[3];
	CSFunctionParser* m_ModeParser[2];

	unsigned int m_numLines[2];
	double** m_ModeDist[2];
};

#endif // PROCESSMODEMATCH_H
