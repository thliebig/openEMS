/*
*	Copyright (C) 2012-2014 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef NF2FF_H
#define NF2FF_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <complex>

using namespace std;

#if defined(WIN32)
	#ifdef BUILD_NF2FF_LIB
	#define NF2FF_EXPORT __declspec(dllexport)
	#else
	#define NF2FF_EXPORT __declspec(dllimport)
	#endif
#else
#define NF2FF_EXPORT
#endif


class TiXmlElement;
class nf2ff_calc;

class NF2FF_EXPORT nf2ff
{
public:
	nf2ff(vector<float> freq, vector<float> theta, vector<float> phi, vector<float> center, unsigned int numThreads=0);
	~nf2ff();

	bool AnalyseFile(string E_Field_file, string H_Field_file);

	void SetRadius(float radius);
	void SetPermittivity(vector<float> permittivity);
	void SetPermeability(vector<float> permeability);
	
	void SetMirror(int type, int dir, float pos);

	double GetTotalRadPower(size_t f_idx) const;
	double GetMaxDirectivity(size_t f_idx) const;

	complex<double>** GetETheta(size_t f_idx) const;
	complex<double>** GetEPhi(size_t f_idx) const;
	double** GetRadPower(size_t f_idx) const;

	//! Write results to a hdf5 file
	bool Write2HDF5(string filename);

	void SetVerboseLevel(int level) {m_Verbose=level;}

	static bool AnalyseXMLNode(TiXmlElement* ti_nf2ff);
	static bool AnalyseXMLFile(string filename);

protected:
	vector<float> m_freq;
	vector<float> m_permittivity;
	vector<float> m_permeability;
	unsigned int m_numTheta;
	unsigned int m_numPhi;
	float* m_theta;
	float* m_phi;
	float m_radius;
	int m_Verbose;
	vector<nf2ff_calc*> m_nf2ff;
};

#endif // NF2FF_H
