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
#include "fparser.hh"
#include "tinyxml.h"
#include "excitation.h"

Excitation::Excitation( double timestep )
{
	Signal_volt = 0;
	Signal_curr = 0;
	E_delay = 0;
	E_amp = 0;
	E_dir = 0;
	Curr_delay = 0;
	Curr_amp = 0;
	Curr_dir = 0;

	for (int n=0;n<3;++n) {
		E_index[n] = 0;
		Curr_index[n] = 0;
	}

	dT = timestep;
	m_nyquistTS = 0;
}

Excitation::~Excitation()
{
	delete[] Signal_volt;
	delete[] Signal_curr;
	delete[] E_delay;
	delete[] E_dir;
	delete[] E_amp;
	delete[] Curr_delay;
	delete[] Curr_dir;
	delete[] Curr_amp;
	for (int n=0;n<3;++n) {
		delete[] E_index[n];
		delete[] Curr_index[n];
	}

}

bool Excitation::setupExcitation( TiXmlElement* Excite, unsigned int maxTS )
{
	if (!Excite) {
		cerr << "Can't read openEMS excitation settings... " << endl;
		return false;
	}

	int Excit_Type=0;
	double f0=0;
	double fc=0;
	Excite->QueryIntAttribute("Type",&Excit_Type);

	switch (Excit_Type)
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
	}

	if (GetNyquistNum() == 0) {
		cerr << "openEMS: excitation setup failed!!" << endl;
		return false;
	}

	return true;
}

unsigned int Excitation::CalcNyquistNum(double fmax) const
{
	if (dT==0) return 1;
	double T0 = 1/fmax;
	return floor(T0/2/dT);
}

unsigned int Excitation::GetMaxExcitationTimestep() const
{
	FDTD_FLOAT maxAmp=0;
	unsigned int maxStep=0;
	for (unsigned int n=1;n<Length+1;++n)
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
	if ((int)Length>nTS)
	{
		cerr << "Operator::CalcGaussianPulsExcitation: Requested excitation pusle would be " << Length << " timesteps or " << Length * dT << " s long. Cutting to max number of timesteps!" << endl;
		Length=(unsigned int)nTS;
	}
//	cerr << "Operator::CalcGaussianPulsExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] Signal_volt;
	delete[] Signal_curr;
	Signal_volt = new FDTD_FLOAT[Length+1];
	Signal_curr = new FDTD_FLOAT[Length+1];
	Signal_volt[0]=0.0;
	Signal_curr[0]=0.0;
	for (unsigned int n=1;n<Length+1;++n)
	{
		double t = (n-1)*dT;
		Signal_volt[n] = cos(2.0*PI*f0*(t-9.0/(2.0*PI*fc)))*exp(-1*pow(2.0*PI*fc*t/3.0-3,2));
		t += 0.5*dT;
		Signal_curr[n] = cos(2.0*PI*f0*(t-9.0/(2.0*PI*fc)))*exp(-1*pow(2.0*PI*fc*t/3.0-3,2));
	}

	SetNyquistNum( CalcNyquistNum(f0+fc) );
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
	for (unsigned int n=1;n<Length+1;++n)
	{
		vars[0] = (n-1)*dT;
		Signal_volt[n] = fParse.Eval(vars);
		vars[0] += 0.5*dT;
		Signal_curr[n] = fParse.Eval(vars);
	}

	SetNyquistNum( CalcNyquistNum(f0) );
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
	for (unsigned int n=1;n<Length+1;++n)
	{
		double t = (n-1)*dT;
		Signal_volt[n] = sin(2.0*PI*f0*t);
		t += 0.5*dT;
		Signal_curr[n] = sin(2.0*PI*f0*t);
	}

	SetNyquistNum( CalcNyquistNum(f0) );
}

void Excitation::setupVoltageExcitation( vector<unsigned int> const volt_vIndex[3], vector<FDTD_FLOAT> const& volt_vExcit,
										 vector<unsigned int> const& volt_vDelay, vector<unsigned int> const& volt_vDir )
{
	E_Count = volt_vIndex[0].size();
	for (int n=0; n<3; n++) {
		delete[] E_index[n];
		E_index[n] = new unsigned int[E_Count];
	}
	delete[] E_delay;
	delete[] E_amp;
	delete[] E_dir;
	E_delay = new unsigned int[E_Count];
	E_amp = new FDTD_FLOAT[E_Count];
	E_dir = new unsigned short[E_Count];

	cerr << "Excitation::setupVoltageExcitation(): Number of voltage excitation points: " << E_Count << endl;
//	if (E_Count==0)
//		cerr << "No E-Field/voltage excitation found!" << endl;
	for (int n=0; n<3; n++)
		for (unsigned int i=0; i<E_Count; i++)
			E_index[n][i] = volt_vIndex[n].at(i);
	for (unsigned int i=0; i<E_Count; i++) {
		E_delay[i] = volt_vDelay.at(i);
		E_amp[i]   = volt_vExcit.at(i);
		E_dir[i]   = volt_vDir.at(i);
	}
}

void Excitation::setupCurrentExcitation( vector<unsigned int> const curr_vIndex[3], vector<FDTD_FLOAT> const& curr_vExcit,
										 vector<unsigned int> const& curr_vDelay, vector<unsigned int> const& curr_vDir )
{
	Curr_Count = curr_vIndex[0].size();
	for (int n=0; n<3; n++) {
		delete[] Curr_index[n];
		Curr_index[n] = new unsigned int[Curr_Count];
	}
	delete[] Curr_delay;
	delete[] Curr_amp;
	delete[] Curr_dir;
	Curr_delay = new unsigned int[Curr_Count];
	Curr_amp = new FDTD_FLOAT[Curr_Count];
	Curr_dir = new unsigned short[Curr_Count];

	cerr << "Excitation::setupCurrentExcitation(): Number of current excitation points: " << Curr_Count << endl;
//	if (Curr_Count==0)
//		cerr << "No H-Field/current excitation found!" << endl;
	for (int n=0;n<3;++n)
		for (unsigned int i=0; i<Curr_Count; i++)
			Curr_index[n][i] = curr_vIndex[n].at(i);
	for (unsigned int i=0; i<Curr_Count; i++) {
		Curr_delay[i] = curr_vDelay.at(i);
		Curr_amp[i]   = curr_vExcit.at(i);
		Curr_dir[i]   = curr_vDir.at(i);
	}

}
