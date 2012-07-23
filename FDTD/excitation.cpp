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

#include "tools/array_ops.h"
#include "tools/useful.h"
#include <iostream>
#include <fstream>
#include "fparser.hh"
#include "tinyxml.h"
#include "excitation.h"

Excitation::Excitation()
{
	Signal_volt = 0;
	Signal_curr = 0;

	m_Excit_Type = -1;

	dT = 0;
	m_nyquistTS = 0;
	m_f_max = 0;
	m_foi = 0;
}

Excitation::~Excitation()
{
	delete[] Signal_volt;
	delete[] Signal_curr;
}

void Excitation::Reset( double timestep )
{
	delete[] Signal_volt;
	Signal_volt = 0;
	delete[] Signal_curr;
	Signal_curr = 0;

	m_Excit_Type = -1;

	dT = timestep;
	m_nyquistTS = 0;
	m_f_max = 0;
	m_foi = 0;
}

bool Excitation::setupExcitation( TiXmlElement* Excite, unsigned int maxTS )
{
	if (!Excite)
	{
		cerr << "Excitation::setupExcitation: Error, can't read openEMS excitation settings... " << endl;
		return false;
	}

	if (dT<=0)
	{
		cerr << "Excitation::setupExcitation: Error, invalid timestep... " << endl;
		return false;
	}

	double f0=0;
	double fc=0;
	Excite->QueryIntAttribute("Type",&m_Excit_Type);

	switch (m_Excit_Type)
	{
	case 0:
		Excite->QueryDoubleAttribute("f0",&f0);
		Excite->QueryDoubleAttribute("fc",&fc);
		CalcGaussianPulsExcitation(f0,fc,maxTS);
		break;
	case 1:
		Excite->QueryDoubleAttribute("f0",&f0);
		CalcSinusExcitation(f0,maxTS);
		break;
	case 2:
		CalcDiracPulsExcitation();
		break;
	case 3:
		CalcStepExcitation();
		break;
	case 10:
		Excite->QueryDoubleAttribute("f0",&f0);
		CalcCustomExcitation(f0,maxTS,Excite->Attribute("Function"));
		break;
	default:
		cerr << "Excitation::setupExcitation: Unknown excitation type: \"" << m_Excit_Type<< "\" !!" << endl;
		m_Excit_Type = -1;
		return false;
	}

	if (GetNyquistNum() == 0)
	{
		cerr << "Excitation::setupExcitation: Unknown error... excitation setup failed!!" << endl;
		return false;
	}

	return true;
}

unsigned int Excitation::GetMaxExcitationTimestep() const
{
	FDTD_FLOAT maxAmp=0;
	unsigned int maxStep=0;
	for (unsigned int n=1; n<Length+1; ++n)
	{
		if (fabs(Signal_volt[n])>maxAmp)
		{
			maxAmp = fabs(Signal_volt[n]);
			maxStep = n;
		}
	}
	return maxStep;
}

void Excitation::CalcGaussianPulsExcitation(double f0, double fc, int nTS)
{
	if (dT==0) return;

	Length = (unsigned int)(2.0 * 9.0/(2.0*PI*fc) / dT);
	if (Length>(unsigned int)nTS)
	{
		cerr << "Operator::CalcGaussianPulsExcitation: Requested excitation pusle would be " << Length << " timesteps or " << Length * dT << " s long. Cutting to max number of timesteps!" << endl;
		Length=(unsigned int)nTS;
	}
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length+1];
	Signal_curr = new FDTD_FLOAT[Length+1];
	Signal_volt[0]=0.0;
	Signal_curr[0]=0.0;
	for (unsigned int n=1; n<Length+1; ++n)
	{
		double t = (n-1)*dT;
		Signal_volt[n] = cos(2.0*PI*f0*(t-9.0/(2.0*PI*fc)))*exp(-1*pow(2.0*PI*fc*t/3.0-3,2));
		t += 0.5*dT;
		Signal_curr[n] = cos(2.0*PI*f0*(t-9.0/(2.0*PI*fc)))*exp(-1*pow(2.0*PI*fc*t/3.0-3,2));
	}

	m_foi = f0;
	m_f_max = f0+fc;

	SetNyquistNum( CalcNyquistNum(f0+fc,dT) );
}

void Excitation::CalcDiracPulsExcitation()
{
	if (dT==0) return;

	Length = 1;
//	cerr << "Operator::CalcDiracPulsExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length+1];
	Signal_curr = new FDTD_FLOAT[Length+1];
	Signal_volt[0]=0.0;
	Signal_volt[1]=1.0;
	Signal_curr[0]=0.0;
	Signal_curr[1]=1.0;

	m_foi = 0;
	m_f_max = 0;

	SetNyquistNum( 1 );
}

void Excitation::CalcStepExcitation()
{
	if (dT==0) return;

	Length = 1;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length+1];
	Signal_curr = new FDTD_FLOAT[Length+1];
	Signal_volt[0]=1.0;
	Signal_volt[1]=1.0;
	Signal_curr[0]=1.0;
	Signal_curr[1]=1.0;

	m_foi = 0;
	m_f_max = 0;

	SetNyquistNum( 1 );
}

void Excitation::CalcCustomExcitation(double f0, int nTS, string signal)
{
	if (dT==0) return;
	if (nTS<=0) return;

	Length = (unsigned int)(nTS);
//	cerr << "Operator::CalcSinusExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length+1];
	Signal_curr = new FDTD_FLOAT[Length+1];
	Signal_volt[0]=0.0;
	Signal_curr[0]=0.0;
	FunctionParser fParse;
	fParse.AddConstant("pi", 3.14159265358979323846);
	fParse.AddConstant("e", 2.71828182845904523536);
	fParse.Parse(signal,"t");
	if (fParse.GetParseErrorType()!=FunctionParser::FP_NO_ERROR)
	{
		cerr << "Operator::CalcCustomExcitation: Function Parser error: " << fParse.ErrorMsg() << endl;
		exit(1);
	}
	double vars[1];
	for (unsigned int n=1; n<Length+1; ++n)
	{
		vars[0] = (n-1)*dT;
		Signal_volt[n] = fParse.Eval(vars);
		vars[0] += 0.5*dT;
		Signal_curr[n] = fParse.Eval(vars);
	}

	m_f_max = f0;
	m_foi = f0;
	SetNyquistNum( CalcNyquistNum(f0,dT) );
}

void Excitation::CalcSinusExcitation(double f0, int nTS)
{
	if (dT==0) return;
	if (nTS<=0) return;

	Length = (unsigned int)(nTS);
//	cerr << "Operator::CalcSinusExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length+1];
	Signal_curr = new FDTD_FLOAT[Length+1];
	Signal_volt[0]=0.0;
	Signal_curr[0]=0.0;
	for (unsigned int n=1; n<Length+1; ++n)
	{
		double t = (n-1)*dT;
		Signal_volt[n] = sin(2.0*PI*f0*t);
		t += 0.5*dT;
		Signal_curr[n] = sin(2.0*PI*f0*t);
	}
	m_f_max = f0;
	m_foi = f0;
	SetNyquistNum( CalcNyquistNum(f0,dT) );
}

void Excitation::DumpVoltageExcite(string filename)
{
	ofstream file;
	file.open( filename.c_str() );
	if (file.fail())
		return;
	for (unsigned int n=1; n<Length+1; ++n)
		file << (n-1)*dT << "\t" << Signal_volt[n] << "\n";
	file.close();
}

void Excitation::DumpCurrentExcite(string filename)
{
	ofstream file;
	file.open( filename.c_str() );
	if (file.fail())
		return;
	for (unsigned int n=1; n<Length+1; ++n)
		file << (n-1)*dT + 0.5*dT << "\t" << Signal_curr[n] << "\n";
	file.close();
}

