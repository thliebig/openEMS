/*
*	Copyright (C) 2012 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef NF2FF_CALC_H
#define NF2FF_CALC_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <complex>
#include <boost/thread.hpp>

using namespace std;

class nf2ff_calc;

// data structure to exchange data between thread-controller and worker-threads
typedef struct
{
	//local working data IN
	int ny;
	float* normDir;
	unsigned int* numLines;
	float **lines;
	float* edge_length_P;
	float* edge_length_PP;

	complex<float>**** E_field;
	complex<float>**** H_field;
	complex<float>**** Js;
	complex<float>**** Ms;

	//local working data OUT
	complex<float>** m_Nt;
	complex<float>** m_Np;
	complex<float>** m_Lt;
	complex<float>** m_Lp;

} nf2ff_data;

class nf2ff_calc_thread
{
public:
	nf2ff_calc_thread(nf2ff_calc* nfc, unsigned int start, unsigned int stop, unsigned int threadID, nf2ff_data &data);
	void operator()();

protected:
	unsigned int m_start, m_stop, m_threadID;
	nf2ff_calc *m_nf_calc;

	nf2ff_data m_data;
};

class nf2ff_calc
{
	// allow full data access to nf2ff_calc_thread class
	friend class nf2ff_calc_thread;
public:
	nf2ff_calc(float freq, vector<float> theta, vector<float> phi);
	~nf2ff_calc();

	float GetRadPower() const {return m_radPower;}
	float GetMaxDirectivity() const {return m_maxDir;}

	complex<float>** GetETheta() const {return m_E_theta;}
	complex<float>** GetEPhi() const {return m_E_phi;}

	unsigned int GetNumThreads() const {return m_numThreads;}
	void SetNumThreads(unsigned int n) {m_numThreads=n;}

	bool AddPlane(float **lines, unsigned int* numLines, complex<float>**** E_field, complex<float>**** H_field);

protected:
	float m_freq;
	float m_radius;

	float m_radPower;
	float m_maxDir;

	complex<float>** m_E_theta;
	complex<float>** m_E_phi;
	complex<float>** m_H_theta;
	complex<float>** m_H_phi;
	float** m_P_rad;

	float m_centerCoord[3];
	unsigned int m_numTheta;
	unsigned int m_numPhi;
	float* m_theta;
	float* m_phi;

	//boost multi-threading
	unsigned int m_numThreads;
	boost::thread_group m_thread_group;
	boost::barrier *m_Barrier;
};


#endif // NF2FF_CALC_H
