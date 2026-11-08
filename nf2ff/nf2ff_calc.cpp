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

#include "nf2ff_calc.h"
#include "../tools/useful.h"
#include "../tools/constants.h"

#include <vector>
#include <cmath>
#include <complex>
#include <iostream>

using namespace std;

// lanes of the complex dot product; the inner loop over them is what the compiler vectorizes
#ifndef NF2FF_LANES
#define NF2FF_LANES 4
#endif

// S[c] += a * sum_q b[q]*f[c*fstride+q] for c<NC, split complex arrays, n a multiple of NF2FF_LANES
template <unsigned int NC>
static inline void nf2ff_cdot(const float* b_re, const float* b_im, const float* f_re, const float* f_im, size_t fstride, unsigned int n, complex<double> a, complex<double>* S)
{
	float s_re[NC][NF2FF_LANES] = {};
	float s_im[NC][NF2FF_LANES] = {};
	for (unsigned int q=0; q<n; q+=NF2FF_LANES)
		for (unsigned int c=0; c<NC; ++c)
			for (unsigned int l=0; l<NF2FF_LANES; ++l)
			{
				float fr = f_re[c*fstride+q+l];
				float fi = f_im[c*fstride+q+l];
				s_re[c][l] += b_re[q+l]*fr - b_im[q+l]*fi;
				s_im[c][l] += b_re[q+l]*fi + b_im[q+l]*fr;
			}
	for (unsigned int c=0; c<NC; ++c)
	{
		double re=0, im=0;
		for (unsigned int l=0; l<NF2FF_LANES; ++l)
		{
			re += s_re[c][l];
			im += s_im[c][l];
		}
		S[c] += a*complex<double>(re, im);
	}
}

// all current components of row p in one pass, so the phase factors are loaded once
static inline void nf2ff_row(const nf2ff_data &d, unsigned int p, const float* b_re, const float* b_im, complex<double> a, complex<double>* S)
{
	size_t off = (size_t)p*d.strideQ;
	size_t fstride = (size_t)d.numP*d.strideQ;
	switch (d.numComp)
	{
	case 4:
		nf2ff_cdot<4>(b_re, b_im, d.F_re+off, d.F_im+off, fstride, d.strideQ, a, S);
		break;
	case 6:
		nf2ff_cdot<6>(b_re, b_im, d.F_re+off, d.F_im+off, fstride, d.strideQ, a, S);
		break;
	default:
		for (unsigned int c=0; c<d.numComp; ++c)
			nf2ff_cdot<1>(b_re, b_im, d.F_re+off+c*fstride, d.F_im+off+c*fstride, 0, d.strideQ, a, S+c);
	}
}

nf2ff_calc_thread::nf2ff_calc_thread(nf2ff_calc* nfc, unsigned int start, unsigned int stop, const nf2ff_data &data)
{
	m_nf_calc = nfc;
	m_start = start;
	m_stop = stop;
	m_data = data;
}

void nf2ff_calc_thread::operator()()
{
	const nf2ff_data &d = m_data;
	const unsigned int numPhi = m_nf_calc->m_numPhi;
	const unsigned int nQ = d.strideQ;

	// phase factors exp(jk*r.r_hat): per point along PP, and per point along P if separable
	vector<float> b_re(nQ, 0), b_im(nQ, 0);
	vector<complex<double>> a(d.numP);

	ArrayLib::ArrayIJ<complex<float>> &Nt = *d.Nt;
	ArrayLib::ArrayIJ<complex<float>> &Np = *d.Np;
	ArrayLib::ArrayIJ<complex<float>> &Lt = *d.Lt;
	ArrayLib::ArrayIJ<complex<float>> &Lp = *d.Lp;

	for (unsigned int idx=m_start; idx<=m_stop; ++idx)
	{
		unsigned int tn = idx / numPhi;
		unsigned int pn = idx % numPhi;
		double sinT = sin((double)m_nf_calc->m_theta[tn]);
		double cosT = cos((double)m_nf_calc->m_theta[tn]);
		double sinP = sin((double)m_nf_calc->m_phi[pn]);
		double cosP = cos((double)m_nf_calc->m_phi[pn]);
		double kdir[3] = {d.k*sinT*cosP, d.k*sinT*sinP, d.k*cosT};

		// sum over the surface of exp(jk*r.r_hat) times each current component
		complex<double> S[6];
		if (d.separable)
		{
			for (unsigned int p=0; p<d.numP; ++p)
			{
				const double* r = d.rP + 3*p;
				a[p] = polar(1.0, r[0]*kdir[0] + r[1]*kdir[1] + r[2]*kdir[2]);
			}
			for (unsigned int q=0; q<d.numPP; ++q)
			{
				const double* r = d.rPP + 3*q;
				double arg = r[0]*kdir[0] + r[1]*kdir[1] + r[2]*kdir[2];
				b_re[q] = cos(arg);
				b_im[q] = sin(arg);
			}
			for (unsigned int p=0; p<d.numP; ++p)
				nf2ff_row(d, p, b_re.data(), b_im.data(), a[p], S);
		}
		else
		{
			for (unsigned int p=0; p<d.numP; ++p)
			{
				for (unsigned int q=0; q<d.numPP; ++q)
				{
					const double* r = d.r + 3*((size_t)p*d.numPP + q);
					double arg = r[0]*kdir[0] + r[1]*kdir[1] + r[2]*kdir[2];
					b_re[q] = cos(arg);
					b_im[q] = sin(arg);
				}
				nf2ff_row(d, p, b_re.data(), b_im.data(), 1.0, S);
			}
		}

		complex<double> J[3], M[3];
		for (unsigned int c=0; c<d.numComp; ++c)
		{
			if (d.comp[c]<3)
				J[d.comp[c]] = S[c];
			else
				M[d.comp[c]-3] = S[c];
		}

		// radiation vectors N and L in spherical components
		Nt(tn, pn) = J[0]*cosT*cosP + J[1]*cosT*sinP - J[2]*sinT;
		Np(tn, pn) = J[1]*cosP - J[0]*sinP;
		Lt(tn, pn) = M[0]*cosT*cosP + M[1]*cosT*sinP - M[2]*sinT;
		Lp(tn, pn) = M[1]*cosP - M[0]*sinP;
	}
}



/***********************************************************************/


nf2ff_calc::nf2ff_calc(float freq, vector<float> theta, vector<float> phi, vector<float> center)
{
	m_freq = freq;
	m_permittivity = 1;
	m_permeability = 1;

	m_numTheta = theta.size();
	m_theta = new float[m_numTheta];
	for (size_t n=0;n<m_numTheta;++n)
		m_theta[n]=theta.at(n);

	m_numPhi = phi.size();
	m_phi = new float[m_numPhi];
	for (size_t n=0;n<m_numPhi;++n)
		m_phi[n]=phi.at(n);

	unsigned int numLines[2] = {m_numTheta, m_numPhi};
	m_E_theta = new ArrayLib::ArrayIJ<std::complex<double>>("E_theta", numLines); // Create2DArray<std::complex<double> >(numLines);
	m_E_phi   = new ArrayLib::ArrayIJ<std::complex<double>>("E_phi"  , numLines); // Create2DArray<std::complex<double> >(numLines);
	m_H_theta = new ArrayLib::ArrayIJ<std::complex<double>>("H_theta", numLines); // Create2DArray<std::complex<double> >(numLines);
	m_H_phi   = new ArrayLib::ArrayIJ<std::complex<double>>("H_phi"  , numLines); // Create2DArray<std::complex<double> >(numLines);
	m_P_rad   = new ArrayLib::ArrayIJ<double>("P_rad", numLines); // Create2DArray<double>(numLines);

	if (center.size()==3)
	{
		m_centerCoord[0]=center.at(0);
		m_centerCoord[1]=center.at(1);
		m_centerCoord[2]=center.at(2);
	}
	else if (center.size()>0)
	{
		cerr << "nf2ff_calc::nf2ff_calc: Warning: Center coordinates error, ignoring!" << endl;
		m_centerCoord[0]=m_centerCoord[1]=m_centerCoord[2]=0.0;
	}
	else
		m_centerCoord[0]=m_centerCoord[1]=m_centerCoord[2]=0.0;

	m_radPower = 0;
	m_maxDir = 0;
	m_radius = 1;

	for (int n=0;n<3;++n)
	{
		m_MirrorType[n] = MIRROR_OFF;
		m_MirrorPos[n]  = 0.0;
	}
	m_numThreads = AvailableThreads();
}

nf2ff_calc::~nf2ff_calc()
{
	delete[] m_phi;
	m_phi = NULL;
	delete[] m_theta;
	m_theta = NULL;

	delete m_E_theta;
	m_E_theta = NULL;
	delete m_E_phi;
	m_E_phi = NULL;
	delete m_H_theta;
	m_H_theta = NULL;
	delete m_H_phi;
	m_H_phi = NULL;
	delete m_P_rad;
	m_P_rad = NULL;
}

int nf2ff_calc::GetNormalDir(unsigned int* numLines)
{
	int ny = -1;
	int nP,nPP;
	for (int n=0;n<3;++n)
	{
		nP = (n+1)%3;
		nPP = (n+2)%3;
		if ((numLines[n]==1) && (numLines[nP]>2) && (numLines[nPP]>2))
			ny=n;
	}
	return ny;
}

void nf2ff_calc::SetMirror(int type, int dir, float pos)
{
	if ((dir<0) || (dir>3))
	{
		cerr << "nf2ff_calc::SetMirror: Error, invalid direction!" << endl;
		return;
	}
	if ((type!=MIRROR_PEC) && (type!=MIRROR_PMC))
	{
		cerr << "nf2ff_calc::SetMirror: Error, invalid type!" << endl;
		return;
	}
	m_MirrorType[dir] = static_cast<MirrorType>(type);
	m_MirrorPos[dir] = pos;
}

bool nf2ff_calc::AddMirrorPlane(int n, float **lines, unsigned int* numLines, ArrayLib::ArrayNIJK<std::complex<float>> &E_field, ArrayLib::ArrayNIJK<std::complex<float>> &H_field, int MeshType)
{
	float E_factor[3] = {1,1,1};
	float H_factor[3] = {1,1,1};

	int nP  = (n+1)%3;
	int nPP = (n+2)%3;

	// mirror in ny direction
	for (unsigned int i=0;i<numLines[n];++i)
		lines[n][i] = 2.0*m_MirrorPos[n] - lines[n][i];
	if (m_MirrorType[n]==MIRROR_PEC)
	{
		H_factor[n]  =-1.0;
		E_factor[nP] =-1.0;
		E_factor[nPP]=-1.0;
	}
	else if (m_MirrorType[n]==MIRROR_PMC)
	{
		E_factor[n]  = -1.0;
		H_factor[nP] = -1.0;
		H_factor[nPP]= -1.0;
	}

	for (int d=0;d<3;++d)
		for (unsigned int i=0;i<numLines[0];++i)
			for (unsigned int j=0;j<numLines[1];++j)
				for (unsigned int k=0;k<numLines[2];++k)
				{
					E_field(d, i, j, k) *= E_factor[d];
					H_field(d, i, j, k) *= H_factor[d];
				}

	return this->AddSinglePlane(lines, numLines, E_field, H_field, MeshType);
}

bool nf2ff_calc::AddPlane(float **lines, unsigned int* numLines, ArrayLib::ArrayNIJK<std::complex<float>> &E_field, ArrayLib::ArrayNIJK<std::complex<float>> &H_field, int MeshType)
{
	this->AddSinglePlane(lines, numLines, E_field, H_field, MeshType);

	for (int n=0;n<3;++n)
	{
		int nP  = (n+1)%3;
		int nPP = (n+2)%3;
		// check if a single mirror plane is on
		if ((m_MirrorType[n]!=MIRROR_OFF) && (m_MirrorType[nP]==MIRROR_OFF) && (m_MirrorType[nPP]==MIRROR_OFF))
		{
			this->AddMirrorPlane(n, lines, numLines, E_field, H_field, MeshType);

			for (unsigned int i=0;i<numLines[n];++i)
				lines[n][i] = 2.0*m_MirrorPos[n] - lines[n][i];

			break;
		}
		//check if two planes are on
		else if ((m_MirrorType[n]==MIRROR_OFF) && (m_MirrorType[nP]!=MIRROR_OFF) && (m_MirrorType[nPP]!=MIRROR_OFF))
		{
			this->AddMirrorPlane(nP, lines, numLines, E_field, H_field, MeshType);
			this->AddMirrorPlane(nPP, lines, numLines, E_field, H_field, MeshType);
			this->AddMirrorPlane(nP, lines, numLines, E_field, H_field, MeshType);

			for (unsigned int i=0;i<numLines[nPP];++i)
				lines[nPP][i] = 2.0*m_MirrorPos[nPP] - lines[nPP][i];

			break;
		}
	}
	// check if all planes are on
	if ((m_MirrorType[0]!=MIRROR_OFF) && (m_MirrorType[1]!=MIRROR_OFF) && (m_MirrorType[2]!=MIRROR_OFF))
	{
		this->AddMirrorPlane(0, lines, numLines, E_field, H_field, MeshType);
		this->AddMirrorPlane(1, lines, numLines, E_field, H_field, MeshType);
		this->AddMirrorPlane(0, lines, numLines, E_field, H_field, MeshType);
		this->AddMirrorPlane(2, lines, numLines, E_field, H_field, MeshType);
		this->AddMirrorPlane(0, lines, numLines, E_field, H_field, MeshType);
		this->AddMirrorPlane(1, lines, numLines, E_field, H_field, MeshType);
		this->AddMirrorPlane(0, lines, numLines, E_field, H_field, MeshType);

		for (unsigned int i=0;i<numLines[2];++i)
			lines[2][i] = 2.0*m_MirrorPos[2] - lines[2][i];

	}

	//cleanup E- & H-Fields
	E_field.Reset();
	H_field.Reset();
	return true;
}

bool nf2ff_calc::AddSinglePlane(float **lines, unsigned int* numLines, ArrayLib::ArrayNIJK<std::complex<float>> &E_field, ArrayLib::ArrayNIJK<std::complex<float>> &H_field, int MeshType)
{
	//find normal direction
	int ny = this->GetNormalDir(numLines);
	if (ny<0)
	{
		cerr << "nf2ff_calc::AddPlane: Error can't determine normal direction..." << endl;
		return false;
	}
	int nP  = (ny+1)%3;
	int nPP = (ny+2)%3;

	ArrayLib::ArrayNIJK<std::complex<float>> Js("Js", numLines);
	ArrayLib::ArrayNIJK<std::complex<float>> Ms("Ms", numLines);

	float normDir[3]= {0,0,0};
	if (lines[ny][0]>=m_centerCoord[ny])
		normDir[ny]=1;
	else
		normDir[ny]=-1;
	unsigned int pos[3];

	float* edge_length_P = new float[numLines[nP]];
	for (unsigned int n=1;n<numLines[nP]-1;++n)
		edge_length_P[n]=0.5*fabs(lines[nP][n+1]-lines[nP][n-1]);
	edge_length_P[0]=0.5*fabs(lines[nP][1]-lines[nP][0]);
	edge_length_P[numLines[nP]-1]=0.5*fabs(lines[nP][numLines[nP]-1]-lines[nP][numLines[nP]-2]);

	float* edge_length_PP = new float[numLines[nPP]];
	for (unsigned int n=1;n<numLines[nPP]-1;++n)
		edge_length_PP[n]=0.5*fabs(lines[nPP][n+1]-lines[nPP][n-1]);
	edge_length_PP[0]=0.5*fabs(lines[nPP][1]-lines[nPP][0]);
	edge_length_PP[numLines[nPP]-1]=0.5*fabs(lines[nPP][numLines[nPP]-1]-lines[nPP][numLines[nPP]-2]);

	//check for cylindrical mesh
	if (MeshType==1)
	{
		if (ny==0) //surface a-z
		{
			for (unsigned int n=0;n<numLines[nP];++n)
				edge_length_P[n]*=lines[0][0]; //angle-width * radius
		}
		else if (ny==2) //surface r-a
		{
			//calculate: area = delta_angle * delta_radius * center_radius
			for (unsigned int n=1;n<numLines[nP]-1;++n)
				edge_length_P[n]*=lines[nP][n];  //radius-width * center-radius
			edge_length_P[0]*=(lines[nP][0]+0.5*edge_length_P[0]);
			edge_length_P[numLines[nP]-1]*=(lines[nP][numLines[nP]-1]-0.5*edge_length_P[numLines[nP]-1]);
		}
	}

	complex<double> power = 0;
	double area;
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				area = edge_length_P[pos[nP]]*edge_length_PP[pos[nPP]];
				power = (E_field(nP, pos[0], pos[1], pos[2])*conj(H_field(nPP, pos[0], pos[1], pos[2])) \
						 - E_field(nPP, pos[0], pos[1], pos[2])*conj(H_field(nP, pos[0], pos[1], pos[2])));
				m_radPower += 0.5*area*real(power)*normDir[ny];
			}
	unsigned int numAngles[2] = {m_numTheta, m_numPhi};

	nf2ff_data data;
	data.numP = numLines[nP];
	data.numPP = numLines[nPP];
	data.strideQ = (data.numPP + NF2FF_LANES-1) / NF2FF_LANES * NF2FF_LANES;
	data.k = 2*PI*m_freq/C0*sqrt(m_permittivity*m_permeability);

	// current components that are not zero by construction: the normal one
	// in Cartesian coordinates, z on a cylindrical r-a surface
	data.numComp = 0;
	for (int c=0; c<6; ++c)
		if ((MeshType==0 && c%3!=ny) || (MeshType==1 && (c%3!=2 || ny!=2)))
			data.comp[data.numComp++] = c;

	// area-weighted surface currents Js = n x H and Ms = -n x E (eq. 8.15a/b),
	// in Cartesian coordinates
	size_t numF = (size_t)data.numComp*data.numP*data.strideQ;
	vector<float> F_re(numF, 0), F_im(numF, 0);
	pos[ny] = 0;
	for (unsigned int p=0; p<data.numP; ++p)
		for (unsigned int q=0; q<data.numPP; ++q)
		{
			pos[nP] = p;
			pos[nPP] = q;
			complex<float> JM[6];
			JM[0] = normDir[1]*H_field(2, pos[0], pos[1], pos[2]) - normDir[2]*H_field(1, pos[0], pos[1], pos[2]);
			JM[1] = normDir[2]*H_field(0, pos[0], pos[1], pos[2]) - normDir[0]*H_field(2, pos[0], pos[1], pos[2]);
			JM[2] = normDir[0]*H_field(1, pos[0], pos[1], pos[2]) - normDir[1]*H_field(0, pos[0], pos[1], pos[2]);
			JM[3] = normDir[2]*E_field(1, pos[0], pos[1], pos[2]) - normDir[1]*E_field(2, pos[0], pos[1], pos[2]);
			JM[4] = normDir[0]*E_field(2, pos[0], pos[1], pos[2]) - normDir[2]*E_field(0, pos[0], pos[1], pos[2]);
			JM[5] = normDir[1]*E_field(0, pos[0], pos[1], pos[2]) - normDir[0]*E_field(1, pos[0], pos[1], pos[2]);
			if (MeshType==1)
			{
				float cos_a = cos(lines[1][pos[1]]);
				float sin_a = sin(lines[1][pos[1]]);
				for (int m=0; m<6; m+=3)
				{
					complex<float> r = JM[m], a = JM[m+1];
					JM[m]   = r*cos_a - a*sin_a;
					JM[m+1] = r*sin_a + a*cos_a;
				}
			}
			float pq_area = edge_length_P[p]*edge_length_PP[q];
			for (unsigned int c=0; c<data.numComp; ++c)
			{
				size_t idx = ((size_t)c*data.numP + p)*data.strideQ + q;
				F_re[idx] = pq_area*real(JM[data.comp[c]]);
				F_im[idx] = pq_area*imag(JM[data.comp[c]]);
			}
		}
	data.F_re = F_re.data();
	data.F_im = F_im.data();

	// point positions relative to the center; separable into a P and a PP part
	// except on a cylindrical r-a surface
	double center[3] = {m_centerCoord[0], m_centerCoord[1], m_centerCoord[2]};
	if (MeshType==1)
	{
		center[0] = m_centerCoord[0]*cos(m_centerCoord[1]);
		center[1] = m_centerCoord[0]*sin(m_centerCoord[1]);
	}
	vector<double> rP, rPP, r;
	data.separable = (MeshType==0) || (ny!=2);
	if (data.separable)
	{
		rP.assign(3*data.numP, 0);
		rPP.assign(3*data.numPP, 0);
		for (unsigned int p=0; p<data.numP; ++p)
		{
			double* v = &rP[3*p];
			if (MeshType==0)
			{
				v[ny]  = lines[ny][0] - center[ny];
				v[nP]  = lines[nP][p] - center[nP];
				v[nPP] = -center[nPP];
			}
			else if (ny==0) // P = alpha, PP = z
			{
				v[0] = lines[0][0]*cos(lines[1][p]) - center[0];
				v[1] = lines[0][0]*sin(lines[1][p]) - center[1];
				v[2] = -center[2];
			}
			else // ny==1: P = z, PP = r
			{
				v[0] = -center[0];
				v[1] = -center[1];
				v[2] = lines[2][p] - center[2];
			}
		}
		for (unsigned int q=0; q<data.numPP; ++q)
		{
			double* v = &rPP[3*q];
			if (MeshType==0)
				v[nPP] = lines[nPP][q];
			else if (ny==0)
				v[2] = lines[2][q];
			else
			{
				v[0] = lines[0][q]*cos(lines[1][0]);
				v[1] = lines[0][q]*sin(lines[1][0]);
			}
		}
	}
	else // cylindrical, ny==2: P = r, PP = alpha
	{
		r.resize(3*(size_t)data.numP*data.numPP);
		for (unsigned int p=0; p<data.numP; ++p)
			for (unsigned int q=0; q<data.numPP; ++q)
			{
				double* v = &r[3*((size_t)p*data.numPP + q)];
				v[0] = lines[0][p]*cos(lines[1][q]) - center[0];
				v[1] = lines[0][p]*sin(lines[1][q]) - center[1];
				v[2] = lines[2][0] - center[2];
			}
	}
	data.rP = rP.data();
	data.rPP = rPP.data();
	data.r = r.data();

	ArrayLib::ArrayIJ<complex<float>> Nt("l_Nt", numAngles);
	ArrayLib::ArrayIJ<complex<float>> Np("l_Np", numAngles);
	ArrayLib::ArrayIJ<complex<float>> Lt("l_Lt", numAngles);
	ArrayLib::ArrayIJ<complex<float>> Lp("l_Lp", numAngles);
	data.Nt = &Nt;
	data.Np = &Np;
	data.Lt = &Lt;
	data.Lp = &Lp;

	// threads calc Nt,Np,Lt and Lp, each for its own range of angles
	vector<unsigned int> jpt = AssignJobs2Threads(m_numTheta*m_numPhi, m_numThreads, true);
	boost::thread_group threads;
	unsigned int start = 0;
	for (size_t n=0; n<jpt.size(); ++n)
	{
		threads.create_thread(nf2ff_calc_thread(this, start, start+jpt.at(n)-1, data));
		start += jpt.at(n);
	}
	threads.join_all();

	ArrayLib::ArrayIJ<std::complex<double>>& E_theta = *m_E_theta;
	ArrayLib::ArrayIJ<std::complex<double>>& E_phi   = *m_E_phi;
	ArrayLib::ArrayIJ<std::complex<double>>& H_theta = *m_H_theta;
	ArrayLib::ArrayIJ<std::complex<double>>& H_phi   = *m_H_phi;
	ArrayLib::ArrayIJ<double>& P_rad = *m_P_rad;

	// calc equations 8.23a/b and 8.24a/b
	float k = data.k;
	complex<float> factor(0,k/4.0/PI/m_radius);
	complex<float> f_exp(0,-1*k*m_radius);
	factor *= exp(f_exp);
	float fZ0 = Z0 * sqrt(m_permeability/m_permittivity);
	complex<float> z_wave = fZ0;
	float P_max = 0;
	for (unsigned int tn=0;tn<m_numTheta;++tn)
		for (unsigned int pn=0;pn<m_numPhi;++pn)
		{
			E_theta(tn, pn) -= factor*(Lp(tn, pn) + z_wave*Nt(tn, pn));
			E_phi(tn, pn)   += factor*(Lt(tn, pn) - z_wave*Np(tn, pn));

			H_theta(tn, pn) += factor*(Np(tn, pn) - Lt(tn, pn)/z_wave);
			H_phi(tn, pn)   -= factor*(Nt(tn, pn) + Lp(tn, pn)/z_wave);

			P_rad(tn, pn) = abs((E_theta(tn, pn)*conj(E_theta(tn, pn))+E_phi(tn, pn)*conj(E_phi(tn, pn))))/(2*fZ0);
			if (P_rad(tn, pn)>P_max)
				P_max = P_rad(tn, pn);
		}

	delete[] edge_length_P; edge_length_P=NULL;
	delete[] edge_length_PP; edge_length_PP=NULL;

	// Maximum directivity (unitless) across sampled angles theta and phi
	// P_max (W/m^2), m_radPower (W), m_radius (m)
	m_maxDir = P_max * (4*PI * m_radius*m_radius / m_radPower);

	return true;
}
