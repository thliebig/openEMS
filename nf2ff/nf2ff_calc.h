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

#include <cstddef>
#include <vector>
#include <cmath>
#include <complex>
#include "../tools/arraylib/array_ij.h"
#include "../tools/arraylib/array_nijk.h"
#include "nf2ff_kernels.h"

enum MirrorType { MIRROR_OFF = 0, MIRROR_PEC = 1, MIRROR_PMC = 2 };

// one surface, prepared for the far-field sums; read-only for the threads
struct nf2ff_data
{
	unsigned int numP, numPP;  // surface points along the two tangential directions
	unsigned int strideQ;      // numPP padded to a multiple of NF2FF_STRIDE_ALIGN
	unsigned int numComp;      // number of non-zero current components
	int comp[6];               // their index: 0..2 = Js x/y/z, 3..5 = Ms x/y/z

	// area-weighted Cartesian surface currents, [comp][p][q], zero padded in q
	const float* F_re;
	const float* F_im;

	// the position r-r_center of point (p,q) is rP[p]+rPP[q] if separable, else r[p][q]
	bool separable;
	const double* rP;   // [p][3]
	const double* rPP;  // [q][3]
	const double* r;    // [p][q][3]

	double k;

	nf2ff_row_func row;   // the kernel for this CPU

	ArrayLib::ArrayIJ<std::complex<float>>* Nt;
	ArrayLib::ArrayIJ<std::complex<float>>* Np;
	ArrayLib::ArrayIJ<std::complex<float>>* Lt;
	ArrayLib::ArrayIJ<std::complex<float>>* Lp;
};

class nf2ff_calc
{
public:
	nf2ff_calc(float freq, std::vector<float> theta, std::vector<float> phi, std::vector<float> center);
	~nf2ff_calc();

	void SetRadius(float radius) {m_radius=radius;}
	void SetPermittivity(float permittivity) {m_permittivity=permittivity;}
	void SetPermeability(float permeability) {m_permeability=permeability;}

	double GetTotalRadPower() const {return m_radPower;}
	double GetMaxDirectivity() const {return m_maxDir;}

	ArrayLib::ArrayIJ<std::complex<double>>* GetETheta() const {return m_E_theta;}
	ArrayLib::ArrayIJ<std::complex<double>>* GetEPhi()   const {return m_E_phi;}
	ArrayLib::ArrayIJ<double>* GetRadPower() const {return m_P_rad;}

	unsigned int GetNumThreads() const {return m_numThreads;}
	void SetNumThreads(unsigned int n) {m_numThreads=n;}

	void SetMirror(int type, int dir, float pos);

	bool AddPlane(float **lines, unsigned int* numLines, ArrayLib::ArrayNIJK<std::complex<float>> &E_field, ArrayLib::ArrayNIJK<std::complex<float>> &H_field, int MeshType=0);

protected:
	float m_freq;
	float m_radius;

	float m_permittivity; //relative electric permittivity
	float m_permeability; //relative magnetic permeability

	double m_radPower;
	double m_maxDir;

	ArrayLib::ArrayIJ<std::complex<double>>* m_E_theta;
	ArrayLib::ArrayIJ<std::complex<double>>* m_E_phi;
	ArrayLib::ArrayIJ<std::complex<double>>* m_H_theta;
	ArrayLib::ArrayIJ<std::complex<double>>* m_H_phi;
	ArrayLib::ArrayIJ<double>* m_P_rad;

	float m_centerCoord[3];
	unsigned int m_numTheta;
	unsigned int m_numPhi;
	float* m_theta;
	float* m_phi;

	//mirror settings
	MirrorType m_MirrorType[3];
	float m_MirrorPos[3];

	int GetNormalDir(unsigned int* numLines);
	bool AddSinglePlane(float **lines, unsigned int* numLines, ArrayLib::ArrayNIJK<std::complex<float>> &E_field, ArrayLib::ArrayNIJK<std::complex<float>> &H_field, int MeshType=0);
	bool AddMirrorPlane(int n, float **lines, unsigned int* numLines, ArrayLib::ArrayNIJK<std::complex<float>> &E_field, ArrayLib::ArrayNIJK<std::complex<float>> &H_field, int MeshType=0);

	//! Nt, Np, Lt and Lp of surface \a d for the flattened (theta,phi) indices [begin, end), on one thread
	void CalcAngles(const nf2ff_data &d, size_t begin, size_t end) const;

	unsigned int m_numThreads;
};


#endif // NF2FF_CALC_H
