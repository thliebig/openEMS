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
#include "excitation.h"

using namespace std;

Excitation::Excitation()
{
	Signal_volt = 0;
	Signal_curr = 0;

	this->Reset(0);

	m_Excit_Type = Excitation::UNDEFINED;
	m_SignalPeriod = 0;
}

Excitation::~Excitation()
{
	this->Reset(0);
}

void Excitation::Reset( double timestep )
{
	delete[] Signal_volt;
	Signal_volt = 0;
	delete[] Signal_curr;
	Signal_curr = 0;

	dT = timestep;
	m_nyquistTS = 0;
	m_f_max = 0;
	m_foi = 0;
}

void Excitation::SetupGaussianPulse(double f0, double fc)
{
	m_Excit_Type = Excitation::GaissianPulse;
	m_f0 = f0;
	m_fc = fc;
	m_f_max = f0+fc;
	m_SignalPeriod = 0;
}

void Excitation::SetupSinusoidal(double f0)
{
	m_Excit_Type = Excitation::Sinusoidal;
	m_f0 = f0;
	m_f_max = f0;
	m_SignalPeriod = 1/f0;
}

void Excitation::SetupDiracPulse(double fmax)
{
	m_Excit_Type = Excitation::DiracPulse;
	m_SignalPeriod = 0;
	m_f_max = fmax;
}

void Excitation::SetupStepExcite(double fmax)
{
	m_Excit_Type = Excitation::Step;
	m_SignalPeriod = 0;
	m_f_max = fmax;
}

void Excitation::SetupCustomExcite(string str, double f0, double fmax)
{
	m_Excit_Type = Excitation::CustomExcite;
	m_CustomExc_Str = str;
	m_f0 = f0;
	m_SignalPeriod = 0;
	m_f_max = fmax;
}

bool Excitation::buildExcitationSignal(unsigned int maxTS)
{
	if (dT<=0)
	{
		cerr << "Excitation::setupExcitation: Error, invalid timestep... " << endl;
		return false;
	}

	switch (m_Excit_Type)
	{
	case Excitation::GaissianPulse:
		CalcGaussianPulsExcitation(m_f0,m_fc,maxTS);
		break;
	case Excitation::Sinusoidal:
		CalcSinusExcitation(m_f0,maxTS);
		break;
	case Excitation::DiracPulse:
		CalcDiracPulsExcitation();
		break;
	case Excitation::Step:
		CalcStepExcitation();
		break;
	case Excitation::CustomExcite:
		CalcCustomExcitation(m_f0,maxTS,m_CustomExc_Str);
		break;
	default:
		cerr << "Excitation::buildExcitationSignal: Unknown excitation type: \"" << m_Excit_Type<< "\" !!" << endl;
		m_Excit_Type = Excitation::UNDEFINED;
		return false;
	}

	if (GetNyquistNum() == 0)
	{
		cerr << "Excitation::buildExcitationSignal: Unknown error... excitation setup failed!!" << endl;
		return false;
	}

	return true;
}

unsigned int Excitation::GetMaxExcitationTimestep() const
{
	FDTD_FLOAT maxAmp=0;
	unsigned int maxStep=0;
	for (unsigned int n=0; n<Length; ++n)
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

	Length = (unsigned int)ceil(2.0 * 9.0/(2.0*PI*fc) / dT);
	if (Length>(unsigned int)nTS)
	{
		cerr << "Operator::CalcGaussianPulsExcitation: Requested excitation pusle would be " << Length << " timesteps or " << Length * dT << " s long. Cutting to max number of timesteps!" << endl;
		Length=(unsigned int)nTS;
	}
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length];
	Signal_curr = new FDTD_FLOAT[Length];
	for (unsigned int n=0; n<Length; ++n)
	{
		double t = n*dT;
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

	Length = 2;
//	cerr << "Operator::CalcDiracPulsExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length];
	Signal_curr = new FDTD_FLOAT[Length];
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

	Length = 2;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length];
	Signal_curr = new FDTD_FLOAT[Length];
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
	Signal_volt = new FDTD_FLOAT[Length];
	Signal_curr = new FDTD_FLOAT[Length];
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
	for (unsigned int n=0; n<Length; ++n)
	{
		vars[0] = n*dT;
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

	Length = (unsigned int)round(2.0/f0/dT);
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length];
	Signal_curr = new FDTD_FLOAT[Length];
	Signal_volt[0]=0.0;
	Signal_curr[0]=0.0;
	for (unsigned int n=1; n<Length; ++n)
	{
		double t = n*dT;
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
	for (unsigned int n=0; n<Length; ++n)
		file << n*dT << "\t" << Signal_volt[n] << "\n";
	file.close();
}

void Excitation::DumpCurrentExcite(string filename)
{
	ofstream file;
	file.open( filename.c_str() );
	if (file.fail())
		return;
	for (unsigned int n=0; n<Length; ++n)
		file << n*dT + 0.5*dT << "\t" << Signal_curr[n] << "\n";
	file.close();
}

