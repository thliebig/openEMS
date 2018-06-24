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

class Excitation
{
public:
	enum ExciteTypes {UNDEFINED=-1, GaissianPulse=0, Sinusoidal=1, DiracPulse=2, Step=3, CustomExcite=10};
	Excitation();
	virtual ~Excitation();

	virtual void Reset( double timestep );

	void SetupGaussianPulse(double f0, double fc);
	void SetupSinusoidal(double f0);
	void SetupDiracPulse(double fmax);
	void SetupStepExcite(double fmax);
	void SetupCustomExcite(std::string str, double f0, double fmax);

	double GetCenterFreq() {return m_f0;}
	double GetCutOffFreq() {return m_fc;}
	double GetMaxFreq() {return m_f_max;}

	bool buildExcitationSignal(unsigned int maxTS);

	//! Get the excitation timestep with the (first) max amplitude
	virtual unsigned int GetMaxExcitationTimestep() const;

	void SetNyquistNum(unsigned int nyquist) {m_nyquistTS=nyquist;}
	unsigned int GetNyquistNum() const {return m_nyquistTS;}

	//! Dump voltage excitation signal to ASCII file
	void DumpVoltageExcite(std::string filename);

	//! Dump current excitation signal to ASCII file
	void DumpCurrentExcite(std::string filename);

	//! Get the used timestep
	double GetTimestep() const {return dT;}

	//! Get the type of excitation
	int GetExciteType() const {return m_Excit_Type;}

	//! Get the length of the excitation signal
	unsigned int GetLength() const {return Length;}

	//! Get the max frequeny excited by this signal
	double GetMaxFrequency() const {return m_f_max;}

	//! Get the frequency of interest
	double GetFrequencyOfInterest() const {return m_foi;}

	//! Get the signal period, 0 if not a periodical signal
	double GetSignalPeriod() const {return m_SignalPeriod;}

	FDTD_FLOAT* GetVoltageSignal() const {return Signal_volt;}
	FDTD_FLOAT* GetCurrentSignal() const {return Signal_curr;}

protected:
	double dT;
	unsigned int m_nyquistTS;
	double m_SignalPeriod;
	ExciteTypes m_Excit_Type;

	//Excitation time-signal
	unsigned int Length;
	FDTD_FLOAT* Signal_volt;
	FDTD_FLOAT* Signal_curr;

	// center frequency
	double m_f0;

	// cutoff-frequency (Gaussian pulse only)
	double m_fc;

	std::string m_CustomExc_Str;

	// max frequency
	double m_f_max;
	// frequency of interest
	double m_foi;

	//! Calculate a custom signal
	virtual void CalcCustomExcitation(double f0, int nTS, std::string signal);
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
