/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef PROCESSFIELDPROBE_H
#define PROCESSFIELDPROBE_H

#include "processintegral.h"

class ProcessFieldProbe : public ProcessIntegral
{
public:
	ProcessFieldProbe(Engine_Interface_Base* eng_if, int type=0);
	virtual ~ProcessFieldProbe();

	virtual std::string GetProcessingName() const;

	virtual std::string GetIntegralName(int row) const;

	//! Set the field type (0 electric field, 1 magnetic field)
	void SetFieldType(int type);

	virtual int GetNumberOfIntegrals() const {return 3;}
	virtual double* CalcMultipleIntegrals();

protected:
	int m_ModeFieldType;
};

#endif // PROCESSFIELDPROBE_H
