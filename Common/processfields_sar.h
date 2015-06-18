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

#ifndef PROCESSFIELDS_SAR_H
#define PROCESSFIELDS_SAR_H

#include "processfields_fd.h"

class ProcessFieldsSAR : public ProcessFieldsFD
{
public:
	ProcessFieldsSAR(Engine_Interface_Base* eng_if);
	virtual ~ProcessFieldsSAR();

	virtual void SetDumpType(DumpType type);

	virtual bool NeedConductivity() const;

	virtual std::string GetProcessingName() const {return "SAR dump";}

	virtual void InitProcess();

	virtual int Process();

	virtual void SetSubSampling(unsigned int subSampleRate, int dir=-1);

	virtual void SetOptResolution(double optRes, int dir=-1);

	//! Set to true for using the conductivity found at the center of a cell, or false for E*J instead
	virtual void SetUseCellConductivity(bool val) {m_UseCellKappa=val;}

	virtual void SetSARAveragingMethod(std::string method) {m_SAR_method=method;}

protected:
	virtual void DumpFDData();

	bool m_UseCellKappa;

	std::string m_SAR_method;

	//! frequency domain electric field storage
	std::vector<std::complex<float>****> m_E_FD_Fields;
	//! frequency domain current density storage
	std::vector<std::complex<float>****> m_J_FD_Fields;
};

#endif // PROCESSFIELDS_SAR_H
