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

#ifndef NF2FF_CALC_H
#define NF2FF_CALC_H

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <complex>
#include <boost/thread.hpp>
#define _USE_MATH_DEFINES

class nf2ff_calc;

#define MIRROR_OFF 0
#define MIRROR_PEC 1
#define MIRROR_PMC 2

// data structure to exchange data between thread-controller and worker-threads
typedef struct
{
	//local working data IN
	int ny;
	int mesh_type;
	float* normDir;
	unsigned int* numLines;
	float **lines;
	float* edge_length_P;
	float* edge_length_PP;

	std::complex<float>**** E_field;
	std::complex<float>**** H_field;
	std::complex<float>**** Js;
	std::complex<float>**** Ms;

	//local working data OUT
	std::complex<double>** m_Nt;
	std::complex<double>** m_Np;
	std::complex<double>** m_Lt;
	std::complex<double>** m_Lp;

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
	nf2ff_calc(float freq, std::vector<float> theta, std::vector<float> phi, std::vector<float> center);
	~nf2ff_calc();

	void SetRadius(float radius) {m_radius=radius;}
	void SetPermittivity(float permittivity) {m_permittivity=permittivity;}
	void SetPermeability(float permeability) {m_permeability=permeability;}

	double GetTotalRadPower() const {return m_radPower;}
	double GetMaxDirectivity() const {return m_maxDir;}

	std::complex<double>** GetETheta() const {return m_E_theta;}
	std::complex<double>** GetEPhi() const {return m_E_phi;}
	double** GetRadPower() const {return m_P_rad;}

	unsigned int GetNumThreads() const {return m_numThreads;}
	void SetNumThreads(unsigned int n) {m_numThreads=n;}

	void SetMirror(int type, int dir, float pos);

	bool AddPlane(float **lines, unsigned int* numLines, std::complex<float>**** E_field, std::complex<float>**** H_field, int MeshType=0);

protected:
	float m_freq;
	float m_radius;

	float m_permittivity; //relative electric permittivity
	float m_permeability; //relative magnetic permeability

	double m_radPower;
	double m_maxDir;

	std::complex<double>** m_E_theta;
	std::complex<double>** m_E_phi;
	std::complex<double>** m_H_theta;
	std::complex<double>** m_H_phi;
	double** m_P_rad;

	float m_centerCoord[3];
	unsigned int m_numTheta;
	unsigned int m_numPhi;
	float* m_theta;
	float* m_phi;

	//mirror settings
	bool m_EnableMirror;
	int m_MirrorType[3];
	float m_MirrorPos[3];

	int GetNormalDir(unsigned int* numLines);
	bool AddSinglePlane(float **lines, unsigned int* numLines, std::complex<float>**** E_field, std::complex<float>**** H_field, int MeshType=0);
	bool AddMirrorPlane(int n, float **lines, unsigned int* numLines, std::complex<float>**** E_field, std::complex<float>**** H_field, int MeshType=0);

	//boost multi-threading
	unsigned int m_numThreads;
	boost::thread_group m_thread_group;
	boost::barrier *m_Barrier;
};


#endif // NF2FF_CALC_H
