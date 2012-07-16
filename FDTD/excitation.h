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

	virtual void Reset( double timestep );

	bool setupExcitation( TiXmlElement* Excite, unsigned int maxTS );

	//! Get the excitation timestep with the (first) max amplitude
	virtual unsigned int GetMaxExcitationTimestep() const;

	void SetNyquistNum(unsigned int nyquist) {m_nyquistTS=nyquist;}
	unsigned int GetNyquistNum() const {return m_nyquistTS;}

	//! Dump voltage excitation signal to ASCII file
	void DumpVoltageExcite(string filename);

	//! Dump current excitation signal to ASCII file
	void DumpCurrentExcite(string filename);

	//! Get the used timestep
	double GetTimestep() const {return dT;}

	//! Get the type of excitation
	int GetExciteType() const {return m_Excit_Type;}

	//! Get the length of the excitation signal
	unsigned int GetLength() const {return Length;}

	FDTD_FLOAT* GetVoltageSignal() const {return Signal_volt;}
	FDTD_FLOAT* GetCurrentSignal() const {return Signal_curr;}

protected:
	double dT;
	unsigned int m_nyquistTS;
	int m_Excit_Type;

	//Excitation time-signal
	unsigned int Length;
	FDTD_FLOAT* Signal_volt;
	FDTD_FLOAT* Signal_curr;

	//! Calculate a custom signal
	virtual void CalcCustomExcitation(double f0, int nTS, string signal);
	//! Calculate an excitation with center of \a f0 and the half bandwidth \a fc
	virtual void CalcGaussianPulsExcitation(double f0, double fc, int nTS);
	//! Calculate a sinusoidal excitation with frequency \a f0 and a duration of \a nTS number of timesteps
	virtual void CalcSinusExcitation(double f0, int nTS);
	//! Calculate a dirac impuls excitation
	virtual void CalcDiracPulsExcitation();
	//! Calculate a step excitation
	virtual void CalcStepExcitation();
};

#endif // EXCITATION_H
