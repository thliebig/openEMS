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

#ifndef EXCITATION_H
#define EXCITATION_H

#include <vector>
#include <string>
#include "tools/constants.h"

class TiXmlElement;

class Excitation
{
public:
	Excitation( double timestep );
	virtual ~Excitation();

	bool setupExcitation( TiXmlElement* Excite, unsigned int maxTS );

	//! Get the excitation timestep with the (first) max amplitude
	virtual unsigned int GetMaxExcitationTimestep() const;
	unsigned int CalcNyquistNum(double fmax) const;

	void setupVoltageExcitation( vector<unsigned int> const volt_vIndex[3], vector<FDTD_FLOAT> const& volt_vExcit,
								 vector<unsigned int> const& volt_vDelay, vector<unsigned int> const& volt_vDir );
	void setupCurrentExcitation( vector<unsigned int> const curr_vIndex[3], vector<FDTD_FLOAT> const& curr_vExcit,
								 vector<unsigned int> const& curr_vDelay, vector<unsigned int> const& curr_vDir );

	void SetNyquistNum(unsigned int nyquist) {m_nyquistTS=nyquist;}
	unsigned int GetNyquistNum() const {return m_nyquistTS;}

	//Excitation time-signal
	unsigned int Length;
	FDTD_FLOAT* Signal_volt;
	FDTD_FLOAT* Signal_curr;

	//E-Field/voltage Excitation
	unsigned int E_Count;
	unsigned int* E_index[3];
	unsigned short* E_dir;
	FDTD_FLOAT* E_amp; //represented as edge-voltages!!
	unsigned int* E_delay;

	//H-Field/current Excitation
	unsigned int Curr_Count;
	unsigned int* Curr_index[3];
	unsigned short* Curr_dir;
	FDTD_FLOAT* Curr_amp; //represented as edge-currents!!
	unsigned int* Curr_delay;

protected:
	double dT;
	unsigned int m_nyquistTS;

	//! Calculate a custom signal
	virtual void CalcCustomExcitation(double f0, int nTS, string signal);
	//! Calculate an excitation with center of f0 and the half bandwidth fc
	virtual void CalcGaussianPulsExcitation(double f0, double fc);
	//! Calculate a sinusoidal excitation with frequency f0 and a duration of nTS number of timesteps
	virtual void CalcSinusExcitation(double f0, int nTS);
	//! Calculate a dirac impuls excitation
	virtual void CalcDiracPulsExcitation();
	//! Calculate a step excitation
	virtual void CalcStepExcitation();
};

#endif // EXCITATION_H
