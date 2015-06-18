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

#include "processfieldprobe.h"

using namespace std;

ProcessFieldProbe::ProcessFieldProbe(Engine_Interface_Base* eng_if, int type) : ProcessIntegral(eng_if)
{
	SetFieldType(type);
}

ProcessFieldProbe::~ProcessFieldProbe()
{

}

string ProcessFieldProbe::GetProcessingName() const
{
	if (m_ModeFieldType==0)
		return "electric field probe";
	if (m_ModeFieldType==1)
		return "magnetic field probe";
	return "unknown field probe";
}

string ProcessFieldProbe::GetIntegralName(int row) const
{
	if (row==0)
	{
		if (m_ModeFieldType==0)
			return "Ex/(V/m)";
		if (m_ModeFieldType==1)
			return "Hx/(A/m)";
	}
	if (row==1)
	{
		if (m_ModeFieldType==0)
			return "Ey/(V/m)";
		if (m_ModeFieldType==1)
			return "Hy/(A/m)";
	}
	if (row==2)
	{
		if (m_ModeFieldType==0)
			return "Ez/(V/m)";
		if (m_ModeFieldType==1)
			return "Hz/(A/m)";
	}
	return "unknown";
}

void ProcessFieldProbe::SetFieldType(int type)
{
	if ((type<0) || (type>1))
	{
		cerr << "ProcessFieldProbe::SetFieldType: Error: unknown field type... skipping" << endl;
		Enabled=false;
	}
	m_ModeFieldType = type;
}

double* ProcessFieldProbe::CalcMultipleIntegrals()
{
	m_Eng_Interface->SetInterpolationType(Engine_Interface_Base::NO_INTERPOLATION);

	switch (m_ModeFieldType)
	{
	case 0:
	default:
		m_Eng_Interface->GetEField(start,m_Results);
		break;
	case 1:
		m_Eng_Interface->GetHField(start,m_Results);
		break;
	}
	return m_Results;
}
