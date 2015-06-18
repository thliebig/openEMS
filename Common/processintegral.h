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

#ifndef PROCESSINTEGRAL_H
#define PROCESSINTEGRAL_H

#include "processing.h"

//! Abstract base class for integral parameter processing
/*!
  \todo Weighting is applied equally to all integral parameter --> todo: weighting for each result individually
  */
class ProcessIntegral : public Processing
{
public:
	virtual ~ProcessIntegral();

	virtual void InitProcess();

	virtual std::string GetProcessingName() const = 0;

	virtual void GetNormalDir(int nd) {m_normDir=nd;}

	//! Flush FD data to disk
	virtual void FlushData();

	//! This method can calculate multiple integral parameter and must be overloaded for each derived class. \sa GetNumberOfIntegrals
	/*!
	  This method will store its integral results internally with a size given by GetNumberOfIntegrals()
	  It will return the result for the CalcIntegral() as default.
	  */
	virtual double* CalcMultipleIntegrals();

	//! Get the name of the integral for the given row. The names will be used in the file header.
	virtual std::string GetIntegralName(int row) const = 0;

	//! Number of calculated results produced by this integral processing. \sa CalcMultipleIntegrals
	virtual int GetNumberOfIntegrals() const {return 1;}

	//! This method should calculate the integral parameter and must be overloaded for each derived class
	virtual double CalcIntegral() {return 0;}

	//! This method will write the TD and FD dump files using CalcIntegral() to calculate the integral parameter
	virtual int Process();

protected:
	ProcessIntegral(Engine_Interface_Base* eng_if);

	void Dump_FD_Data(double factor, std::string filename);

	std::vector<double_complex> *m_FD_Results;
	double *m_Results;

	int m_normDir; // normal direction as required by some integral processings
};

#endif // PROCESSINTEGRAL_H
