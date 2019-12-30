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
#include "../tools/array_ops.h"
#include "../tools/useful.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>

using namespace std;

nf2ff_calc_thread::nf2ff_calc_thread(nf2ff_calc* nfc, unsigned int start, unsigned int stop, unsigned int threadID, nf2ff_data &data)
{
	m_nf_calc = nfc;
	m_start = start;
	m_stop = stop;
	m_threadID = threadID;
	m_data = data;
}

void nf2ff_calc_thread::operator()()
{
	m_nf_calc->m_Barrier->wait(); // start

	int ny = m_data.ny;
	int nP = (ny+1)%3;
	int nPP = (ny+2)%3;

	unsigned int* numLines = m_data.numLines;
	float* normDir = m_data.normDir;
	float **lines = m_data.lines;
	float* edge_length_P = m_data.edge_length_P;
	float* edge_length_PP = m_data.edge_length_PP;

	unsigned int pos[3];
	unsigned int pos_t=0;
	unsigned int num_t=m_stop-m_start+1;


	complex<float>**** Js=m_data.Js;
	complex<float>**** Ms=m_data.Ms;
	complex<float>**** E_field=m_data.E_field;
	complex<float>**** H_field=m_data.H_field;

	int mesh_type = m_data.mesh_type;

	// calc Js and Ms (eq. 8.15a/b)
	pos[ny]=0;
	for (pos_t=0; pos_t<num_t; ++pos_t)
	{
		pos[nP] = m_start+pos_t;
		for (pos[nPP]=0; pos[nPP]<numLines[nPP]; ++pos[nPP])
		{
			// Js =  n x H
			Js[0][pos[0]][pos[1]][pos[2]] = normDir[1]*H_field[2][pos[0]][pos[1]][pos[2]] - normDir[2]*H_field[1][pos[0]][pos[1]][pos[2]];
			Js[1][pos[0]][pos[1]][pos[2]] = normDir[2]*H_field[0][pos[0]][pos[1]][pos[2]] - normDir[0]*H_field[2][pos[0]][pos[1]][pos[2]];
			Js[2][pos[0]][pos[1]][pos[2]] = normDir[0]*H_field[1][pos[0]][pos[1]][pos[2]] - normDir[1]*H_field[0][pos[0]][pos[1]][pos[2]];

			// Ms = -n x E
			Ms[0][pos[0]][pos[1]][pos[2]] = normDir[2]*E_field[1][pos[0]][pos[1]][pos[2]] - normDir[1]*E_field[2][pos[0]][pos[1]][pos[2]];
			Ms[1][pos[0]][pos[1]][pos[2]] = normDir[0]*E_field[2][pos[0]][pos[1]][pos[2]] - normDir[2]*E_field[0][pos[0]][pos[1]][pos[2]];
			Ms[2][pos[0]][pos[1]][pos[2]] = normDir[1]*E_field[0][pos[0]][pos[1]][pos[2]] - normDir[0]*E_field[1][pos[0]][pos[1]][pos[2]];

			//transform to cartesian coordinates
			if (mesh_type==1)
			{
				Js[0][pos[0]][pos[1]][pos[2]] = (normDir[1]*H_field[2][pos[0]][pos[1]][pos[2]] - normDir[2]*H_field[1][pos[0]][pos[1]][pos[2]])*cos(lines[1][pos[1]]) \
						- (normDir[2]*H_field[0][pos[0]][pos[1]][pos[2]] - normDir[0]*H_field[2][pos[0]][pos[1]][pos[2]])*sin(lines[1][pos[1]]);
				Js[1][pos[0]][pos[1]][pos[2]] = (normDir[1]*H_field[2][pos[0]][pos[1]][pos[2]] - normDir[2]*H_field[1][pos[0]][pos[1]][pos[2]])*sin(lines[1][pos[1]]) \
						+ (normDir[2]*H_field[0][pos[0]][pos[1]][pos[2]] - normDir[0]*H_field[2][pos[0]][pos[1]][pos[2]])*cos(lines[1][pos[1]]);

				Ms[0][pos[0]][pos[1]][pos[2]] = (normDir[2]*E_field[1][pos[0]][pos[1]][pos[2]] - normDir[1]*E_field[2][pos[0]][pos[1]][pos[2]])*cos(lines[1][pos[1]]) \
						- (normDir[0]*E_field[2][pos[0]][pos[1]][pos[2]] - normDir[2]*E_field[0][pos[0]][pos[1]][pos[2]])*sin(lines[1][pos[1]]);
				Ms[1][pos[0]][pos[1]][pos[2]] = (normDir[2]*E_field[1][pos[0]][pos[1]][pos[2]] - normDir[1]*E_field[2][pos[0]][pos[1]][pos[2]])*sin(lines[1][pos[1]]) \
						+ (normDir[0]*E_field[2][pos[0]][pos[1]][pos[2]] - normDir[2]*E_field[0][pos[0]][pos[1]][pos[2]])*cos(lines[1][pos[1]]);
			}
		}
	}

	complex<double>** m_Nt=m_data.m_Nt;
	complex<double>** m_Np=m_data.m_Np;
	complex<double>** m_Lt=m_data.m_Lt;
	complex<double>** m_Lp=m_data.m_Lp;

	float center[3] = {m_nf_calc->m_centerCoord[0],m_nf_calc->m_centerCoord[1],m_nf_calc->m_centerCoord[2]};
	if (mesh_type==1)
	{
		center[0] = m_nf_calc->m_centerCoord[0]*cos(m_nf_calc->m_centerCoord[1]);
		center[1] = m_nf_calc->m_centerCoord[0]*sin(m_nf_calc->m_centerCoord[1]);
	}
	// calc local Nt,Np,Lt and Lp
	float area;
	float cosT_cosP,cosP_sinT;
	float cosT_sinP,sinT_sinP;
	float sinT,sinP;
	float cosP,cosT;
	float r_cos_psi;
	float k = 2*M_PI*m_nf_calc->m_freq/__C0__*sqrt(m_nf_calc->m_permittivity*m_nf_calc->m_permeability);
	complex<float> exp_jkr;
	complex<float> _I_(0,1);
	for (unsigned int tn=0;tn<m_nf_calc->m_numTheta;++tn)
		for (unsigned int pn=0;pn<m_nf_calc->m_numPhi;++pn)
		{
			sinT = sin(m_nf_calc->m_theta[tn]);
			sinP = sin(m_nf_calc->m_phi[pn]);
			cosT = cos(m_nf_calc->m_theta[tn]);
			cosP = cos(m_nf_calc->m_phi[pn]);
			cosT_cosP = cosT*cosP;
			cosT_sinP = cosT*sinP;
			cosP_sinT = cosP*sinT;
			sinT_sinP = sinP*sinT;

			for (pos_t=0; pos_t<num_t; ++pos_t)
			{
				pos[nP] = m_start+pos_t;
				for (pos[nPP]=0; pos[nPP]<numLines[nPP]; ++pos[nPP])
				{
					if (mesh_type==0)
						r_cos_psi = (lines[0][pos[0]]-center[0])*cosP_sinT + (lines[1][pos[1]]-center[1])*sinT_sinP + (lines[2][pos[2]]-center[2])*cosT;
					else
						r_cos_psi = ((lines[0][pos[0]]*cos(lines[1][pos[1]]))-center[0])*cosP_sinT + ((lines[0][pos[0]]*sin(lines[1][pos[1]]))-center[1])*sinT_sinP + (lines[2][pos[2]]-center[2])*cosT;
					exp_jkr = exp(_I_*k*r_cos_psi);
					area = edge_length_P[pos[nP]]*edge_length_PP[pos[nPP]];
					m_Nt[tn][pn] += area*exp_jkr*(Js[0][pos[0]][pos[1]][pos[2]]*cosT_cosP + Js[1][pos[0]][pos[1]][pos[2]]*cosT_sinP \
												  - Js[2][pos[0]][pos[1]][pos[2]]*sinT);
					m_Np[tn][pn] += area*exp_jkr*(Js[1][pos[0]][pos[1]][pos[2]]*cosP - Js[0][pos[0]][pos[1]][pos[2]]*sinP);

					m_Lt[tn][pn] += area*exp_jkr*(Ms[0][pos[0]][pos[1]][pos[2]]*cosT_cosP + Ms[1][pos[0]][pos[1]][pos[2]]*cosT_sinP \
												  - Ms[2][pos[0]][pos[1]][pos[2]]*sinT);
					m_Lp[tn][pn] += area*exp_jkr*(Ms[1][pos[0]][pos[1]][pos[2]]*cosP - Ms[0][pos[0]][pos[1]][pos[2]]*sinP);
				}
			}
		}

	m_nf_calc->m_Barrier->wait(); //combine all thread local Nt,Np,Lt and Lp

	m_nf_calc->m_Barrier->wait(); //wait for termination
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
	m_E_theta = Create2DArray<std::complex<double> >(numLines);
	m_E_phi = Create2DArray<std::complex<double> >(numLines);
	m_H_theta = Create2DArray<std::complex<double> >(numLines);
	m_H_phi = Create2DArray<std::complex<double> >(numLines);
	m_P_rad = Create2DArray<double>(numLines);

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

	m_Barrier = NULL;
	m_numThreads = boost::thread::hardware_concurrency();
}

nf2ff_calc::~nf2ff_calc()
{
	delete[] m_phi;
	m_phi = NULL;
	delete[] m_theta;
	m_theta = NULL;

	unsigned int numLines[2] = {m_numTheta, m_numPhi};
	Delete2DArray(m_E_theta,numLines);
	m_E_theta = NULL;
	Delete2DArray(m_E_phi,numLines);
	m_E_phi = NULL;
	Delete2DArray(m_H_theta,numLines);
	m_H_theta = NULL;
	Delete2DArray(m_H_phi,numLines);
	m_H_phi = NULL;
	Delete2DArray(m_P_rad,numLines);
	m_P_rad = NULL;

	delete m_Barrier;
	m_Barrier = NULL;
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
	m_MirrorType[dir] = type;
	m_MirrorPos[dir] = pos;
}

bool nf2ff_calc::AddMirrorPlane(int n, float **lines, unsigned int* numLines, complex<float>**** E_field, complex<float>**** H_field, int MeshType)
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
					E_field[d][i][j][k] *= E_factor[d];
					H_field[d][i][j][k] *= H_factor[d];
				}

	return this->AddSinglePlane(lines, numLines, E_field, H_field, MeshType);
}

bool nf2ff_calc::AddPlane(float **lines, unsigned int* numLines, complex<float>**** E_field, complex<float>**** H_field, int MeshType)
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
	Delete_N_3DArray(E_field,numLines);
	Delete_N_3DArray(H_field,numLines);
	return true;
}

bool nf2ff_calc::AddSinglePlane(float **lines, unsigned int* numLines, complex<float>**** E_field, complex<float>**** H_field, int MeshType)
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

	complex<float>**** Js = Create_N_3DArray<complex<float> >(numLines);
	complex<float>**** Ms = Create_N_3DArray<complex<float> >(numLines);

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
				power = (E_field[nP][pos[0]][pos[1]][pos[2]]*conj(H_field[nPP][pos[0]][pos[1]][pos[2]]) \
						 - E_field[nPP][pos[0]][pos[1]][pos[2]]*conj(H_field[nP][pos[0]][pos[1]][pos[2]]));
				m_radPower += 0.5*area*real(power)*normDir[ny];
			}
	unsigned int numAngles[2] = {m_numTheta, m_numPhi};

	// setup multi-threading jobs
	vector<unsigned int> jpt = AssignJobs2Threads(numLines[nP], m_numThreads, true);
	m_numThreads = jpt.size();
	nf2ff_data* thread_data = new nf2ff_data[m_numThreads];
	m_Barrier = new boost::barrier(m_numThreads+1); // numThread workers + 1 controller
	unsigned int start=0;
	unsigned int stop=jpt.at(0)-1;
	for (unsigned int n=0; n<m_numThreads; n++)
	{
		thread_data[n].ny=ny;
		thread_data[n].mesh_type = MeshType;
		thread_data[n].normDir=normDir;
		thread_data[n].numLines=numLines;
		thread_data[n].lines=lines;
		thread_data[n].edge_length_P=edge_length_P;
		thread_data[n].edge_length_PP=edge_length_PP;
		thread_data[n].E_field=E_field;
		thread_data[n].H_field=H_field;
		thread_data[n].Js=Js;
		thread_data[n].Ms=Ms;
		thread_data[n].m_Nt=Create2DArray<complex<double> >(numAngles);
		thread_data[n].m_Np=Create2DArray<complex<double> >(numAngles);
		thread_data[n].m_Lt=Create2DArray<complex<double> >(numAngles);
		thread_data[n].m_Lp=Create2DArray<complex<double> >(numAngles);

		boost::thread *t = new boost::thread( nf2ff_calc_thread(this,start,stop,n,thread_data[n]) );

		m_thread_group.add_thread( t );

		start = stop+1;
		if (n<m_numThreads-1)
			stop = start + jpt.at(n+1)-1;
	}
	//all threads a running and waiting for the barrier

	m_Barrier->wait(); //start

	// threads: calc Js and Ms (eq. 8.15a/b)
	// threads calc their local Nt,Np,Lt and Lp

	m_Barrier->wait(); //combine all thread local Nt,Np,Lt and Lp

	complex<float>** Nt = Create2DArray<complex<float> >(numAngles);
	complex<float>** Np = Create2DArray<complex<float> >(numAngles);
	complex<float>** Lt = Create2DArray<complex<float> >(numAngles);
	complex<float>** Lp = Create2DArray<complex<float> >(numAngles);

	for (unsigned int n=0; n<m_numThreads; n++)
	{
		for (unsigned int tn=0;tn<m_numTheta;++tn)
			for (unsigned int pn=0;pn<m_numPhi;++pn)
			{
				Nt[tn][pn] += thread_data[n].m_Nt[tn][pn];
				Np[tn][pn] += thread_data[n].m_Np[tn][pn];
				Lt[tn][pn] += thread_data[n].m_Lt[tn][pn];
				Lp[tn][pn] += thread_data[n].m_Lp[tn][pn];
			}
		Delete2DArray(thread_data[n].m_Nt,numAngles);
		Delete2DArray(thread_data[n].m_Np,numAngles);
		Delete2DArray(thread_data[n].m_Lt,numAngles);
		Delete2DArray(thread_data[n].m_Lp,numAngles);
	}

	m_Barrier->wait(); //wait for termination
	m_thread_group.join_all(); // wait for termination
	delete m_Barrier;
	m_Barrier = NULL;

	//cleanup Js & Ms
	Delete_N_3DArray(Js,numLines);
	Delete_N_3DArray(Ms,numLines);

	// calc equations 8.23a/b and 8.24a/b
	float k = 2*M_PI*m_freq/__C0__*sqrt(m_permittivity*m_permeability);
	complex<float> factor(0,k/4.0/M_PI/m_radius);
	complex<float> f_exp(0,-1*k*m_radius);
	factor *= exp(f_exp);
	float fZ0 = __Z0__ * sqrt(m_permeability/m_permittivity);
	complex<float> Z0 = fZ0;
	float P_max = 0;
	for (unsigned int tn=0;tn<m_numTheta;++tn)
		for (unsigned int pn=0;pn<m_numPhi;++pn)
		{
			m_E_theta[tn][pn] -= factor*(Lp[tn][pn] + Z0*Nt[tn][pn]);
			m_E_phi[tn][pn] += factor*(Lt[tn][pn] - Z0*Np[tn][pn]);

			m_H_theta[tn][pn] += factor*(Np[tn][pn] - Lt[tn][pn]/Z0);
			m_H_phi[tn][pn] -= factor*(Nt[tn][pn] + Lp[tn][pn]/Z0);

			m_P_rad[tn][pn] = abs((m_E_theta[tn][pn]*conj(m_E_theta[tn][pn])+m_E_phi[tn][pn]*conj(m_E_phi[tn][pn])))/(2*fZ0);
			if (m_P_rad[tn][pn]>P_max)
				P_max = m_P_rad[tn][pn];
		}

	//cleanup Nx and Lx
	Delete2DArray(Nt,numAngles);
	Delete2DArray(Np,numAngles);
	Delete2DArray(Lt,numAngles);
	Delete2DArray(Lp,numAngles);
	delete[] edge_length_P; edge_length_P=NULL;
	delete[] edge_length_PP; edge_length_PP=NULL;
	delete[] thread_data; thread_data=NULL;

	m_maxDir = 4*M_PI*P_max / m_radPower;

	return true;
}
