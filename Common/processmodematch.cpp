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

#include "processmodematch.h"
#include "CSFunctionParser.h"
#include "Common/operator_base.h"
#include "tools/array_ops.h"

using namespace std;

ProcessModeMatch::ProcessModeMatch(Engine_Interface_Base* eng_if) : ProcessIntegral(eng_if)
{
	for (int n=0; n<2; ++n)
	{
		m_ModeParser[n] = new CSFunctionParser();
		m_ModeDist[n] = NULL;
	}
	delete[] m_Results;
	m_Results = new double[2];
}

ProcessModeMatch::~ProcessModeMatch()
{
	for (int n=0; n<2; ++n)
	{
		delete m_ModeParser[n];
		m_ModeParser[n] = NULL;
	}
	Reset();
}

string ProcessModeMatch::GetIntegralName(int row) const
{
	if (row==0)
	{
		if (m_ModeFieldType==0)
			return "voltage";
		if (m_ModeFieldType==1)
			return "current";
	}
	if (row==1)
	{
		return "mode_purity";
	}
	return "unknown";
}

string ProcessModeMatch::GetProcessingName() const
{
	if (m_ModeFieldType==0)
		return "voltage mode matching";
	if (m_ModeFieldType==1)
		return "current mode matching";
	return "unknown mode matching";
}

void ProcessModeMatch::InitProcess()
{
	if (!Enabled) return;

	if (m_Eng_Interface==NULL)
	{
		cerr << "ProcessModeMatch::InitProcess: Error, Engine_Interface is NULL, abort mode mathcing..." << endl;
		Enabled=false;
		return;
	}
	m_Eng_Interface->SetInterpolationType(Engine_Interface_Base::NODE_INTERPOLATE);

	int Dump_Dim=0;
	m_ny = -1;
	for (int n=0; n<3; ++n)
	{
		if (start[n]>stop[n])
		{
			unsigned int help=start[n];
			start[n]=stop[n];
			stop[n]=help;
		}

		//exclude boundaries from mode-matching
		if (start[n]==0)
			++start[n];
		if (stop[n]==Op->GetNumberOfLines(n)-1)
			--stop[n];

		if (stop[n]!=start[n])
			++Dump_Dim;

		if (stop[n] == start[n])
			m_ny = n;
	}

	if ((Dump_Dim!=2) || (m_ny<0))
	{
		cerr << "ProcessModeMatch::InitProcess(): Warning Mode Matching Integration Box \"" << m_filename << "\" is not a surface (found dimension: " << Dump_Dim << ")" << endl;
		SetEnable(false);
		Reset();
		return;
	}

	int nP = (m_ny+1)%3;
	int nPP = (m_ny+2)%3;
	m_numLines[0] = stop[nP] - start[nP] + 1;
	m_numLines[1] = stop[nPP] - start[nPP] + 1;

	for (int n=0; n<2; ++n)
	{
		int ny = (m_ny+n+1)%3;
		int res = m_ModeParser[n]->Parse(m_ModeFunction[ny], "x,y,z,rho,a,r,t");
		if (res >= 0)
		{
			cerr << "ProcessModeMatch::InitProcess(): Warning, an error occured parsing the mode matching function (see below) ..." << endl;
			cerr << m_ModeFunction[ny] << "\n" << string(res, ' ') << "^\n" << m_ModeParser[n]->ErrorMsg() << "\n";
			SetEnable(false);
			Reset();
		}
	}

	for (int n=0; n<2; ++n)
	{
		m_ModeDist[n] = Create2DArray<double>(m_numLines);
	}

	bool dualMesh = m_ModeFieldType==1;
	unsigned int pos[3] = {0,0,0};
	double discLine[3] = {0,0,0};
	double gridDelta = 1; // 1 -> mode-matching function is definied in drawing units...
	double var[7];
	pos[m_ny] = start[m_ny];
	discLine[m_ny] = Op->GetDiscLine(m_ny,pos[m_ny],dualMesh);
	double norm = 0;
	double area = 0;
	for (unsigned int posP = 0; posP<m_numLines[0]; ++posP)
	{
		pos[nP] = start[nP] + posP;
		discLine[nP] = Op->GetDiscLine(nP,pos[nP],dualMesh);
		for (unsigned int posPP = 0; posPP<m_numLines[1]; ++posPP)
		{
			pos[nPP] = start[nPP] + posPP;
			discLine[nPP] = Op->GetDiscLine(nPP,pos[nPP],dualMesh);

			var[0] = discLine[0] * gridDelta; // x
			var[1] = discLine[1] * gridDelta; // y
			var[2] = discLine[2] * gridDelta; // z
			var[3] = sqrt(discLine[0]*discLine[0] + discLine[1]*discLine[1]) * gridDelta; // rho = sqrt(x^2 + y^2)
			var[4] = atan2(discLine[1], discLine[0]); // a = atan(y,x)
			var[5] = sqrt(pow(discLine[0],2)+pow(discLine[1],2)+pow(discLine[2],2)) * gridDelta; // r
			var[6] = asin(1)-atan(var[2]/var[3]); //theta (t)

			if (m_Mesh_Type == CYLINDRICAL_MESH)
			{
				var[3] = discLine[0] * gridDelta; // rho
				var[4] = discLine[1]; // a
				var[0] = discLine[0] * cos(discLine[1]) * gridDelta; // x = r*cos(a)
				var[1] = discLine[0] * sin(discLine[1]) * gridDelta; // y = r*sin(a)
				var[5] = sqrt(pow(discLine[0],2)+pow(discLine[2],2)) * gridDelta; // r
				var[6] = asin(1)-atan(var[2]/var[3]); //theta (t)
			}
			area = Op->GetNodeArea(m_ny,pos,dualMesh);
			for (int n=0; n<2; ++n)
			{
				m_ModeDist[n][posP][posPP] = m_ModeParser[n]->Eval(var); //calc mode template
				if ((std::isnan(m_ModeDist[n][posP][posPP])) || (std::isinf(m_ModeDist[n][posP][posPP])))
					m_ModeDist[n][posP][posPP] = 0.0;
				norm += pow(m_ModeDist[n][posP][posPP],2) * area;
			}
//			cerr << discLine[0] << " " << discLine[1] << " : " << m_ModeDist[0][posP][posPP] << " , " << m_ModeDist[1][posP][posPP] << endl;
		}
	}

	norm = sqrt(norm);
//	cerr << norm << endl;
	// normalize template function...
	for (unsigned int posP = 0; posP<m_numLines[0]; ++posP)
		for (unsigned int posPP = 0; posPP<m_numLines[1]; ++posPP)
		{
			for (int n=0; n<2; ++n)
			{
				m_ModeDist[n][posP][posPP] /= norm;
			}
//			cerr << posP << " " << posPP << " : " << m_ModeDist[0][posP][posPP] << " , " << m_ModeDist[1][posP][posPP] << endl;
		}

	ProcessIntegral::InitProcess();
}

void ProcessModeMatch::Reset()
{
	ProcessIntegral::Reset();
	for (int n=0; n<2; ++n)
	{
		Delete2DArray<double>(m_ModeDist[n],m_numLines);
		m_ModeDist[n] = NULL;
	}
}


void ProcessModeMatch::SetModeFunction(int ny, string function)
{
	if ((ny<0) || (ny>2)) return;
	m_ModeFunction[ny] = function;
}

void ProcessModeMatch::SetFieldType(int type)
{
	m_ModeFieldType = type;
	if ((type<0) || (type>1))
		cerr << "ProcessModeMatch::SetFieldType: Warning, unknown field type..." << endl;
}

double* ProcessModeMatch::CalcMultipleIntegrals()
{
	double value = 0;
	double field = 0;
	double purity = 0;
	double area = 0;
    bool dualMesh = m_ModeFieldType==1;

	int nP = (m_ny+1)%3;
	int nPP = (m_ny+2)%3;

	unsigned int pos[3] = {0,0,0};
	pos[m_ny] = start[m_ny];

	double out[3]={0,0,0};

	for (unsigned int posP = 0; posP<m_numLines[0]; ++posP)
	{
		pos[nP] = start[nP] + posP;
		for (unsigned int posPP = 0; posPP<m_numLines[1]; ++posPP)
		{
			pos[nPP] = start[nPP] + posPP;
			area = Op->GetNodeArea(m_ny,pos,dualMesh);
			if (m_ModeFieldType==0)
				m_Eng_Interface->GetEField(pos,out);
			if (m_ModeFieldType==1)
				m_Eng_Interface->GetHField(pos,out);

			for (int n=0; n<2; ++n)
			{
				field = out[(m_ny+n+1)%3];
				value += field * m_ModeDist[n][posP][posPP] * area;
				purity += field*field * area;
			}
		}
	}
	if (purity!=0)
		m_Results[1] = value*value/purity;
	else
		m_Results[1] = 0;
	m_Results[0] = value;
	return m_Results;
}
