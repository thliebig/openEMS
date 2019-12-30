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


#include <algorithm>
#include "sar_calculation.h"
#include "cfloat"
#include "array_ops.h"
#include "global.h"

using namespace std;

SAR_Calculation::SAR_Calculation()
{
	m_Vx_Used = NULL;
	m_Vx_Valid = NULL;
	m_DebugLevel = 0;
	SetAveragingMethod(SIMPLE, true);
	Reset();
}

void SAR_Calculation::Reset()
{
	Delete3DArray(m_Vx_Used,m_numLines);
	m_Vx_Used = NULL;
	Delete3DArray(m_Vx_Valid,m_numLines);
	m_Vx_Valid = NULL;

	m_avg_mass = 0;
	m_numLines[0]=m_numLines[1]=m_numLines[2]=0;
	m_cellWidth[0]=m_cellWidth[1]=m_cellWidth[2]=NULL;

	m_cell_volume = NULL;
	m_cell_density = NULL;
	m_cell_conductivity = NULL;
	m_E_field = NULL;
	m_J_field = NULL;

	Delete3DArray(m_Vx_Used,m_numLines);
	m_Vx_Used = NULL;
	Delete3DArray(m_Vx_Valid,m_numLines);
	m_Vx_Valid = NULL;
}

void SAR_Calculation::SetAveragingMethod(string method, bool silent)
{
	if (method.compare("IEEE_C95_3")==0)
		return SetAveragingMethod(IEEE_C95_3, silent);
	if (method.compare("IEEE_62704")==0)
		return SetAveragingMethod(IEEE_62704, silent);
	if (method.compare("Simple")==0)
		return SetAveragingMethod(SIMPLE, silent);

	cerr << __func__ << ": Error, """ << method << """ is an unknown averaging method..." << endl;
	// unknown method, fallback to simple
	SetAveragingMethod(SIMPLE, false);
}

void SAR_Calculation::SetAveragingMethod(SARAveragingMethod method, bool silent)
{
	if (method==IEEE_62704)
	{
		m_massTolerance = 0.0001;
		m_maxMassIterations = 100;
		m_maxBGRatio = 0.1;
		m_markPartialAsUsed = false;
		m_UnusedRelativeVolLimit = 1.05;
		m_IgnoreFaceValid = false;
		if (!silent)
			cerr << __func__ << ": Setting averaging method to IEEE_62704" << endl;
		return;
	}
	else if (method==IEEE_C95_3)
	{
		m_massTolerance = 0.05;
		m_maxMassIterations = 100;
		m_maxBGRatio = 1;
		m_markPartialAsUsed = true;
		m_UnusedRelativeVolLimit = 1;
		m_IgnoreFaceValid = false;
		if (!silent)
			cerr << __func__ << ": Setting averaging method to IEEE_C95_3" << endl;
		return;
	}
	else if (method==SIMPLE)
	{
		m_massTolerance = 0.05;
		m_maxMassIterations = 100;
		m_maxBGRatio = 1;
		m_markPartialAsUsed = true;
		m_UnusedRelativeVolLimit = 1;
		m_IgnoreFaceValid = true;
		if (!silent)
			cerr << __func__ << ": Setting averaging method to SIMPLE" << endl;
		return;
	}

	cerr << __func__ << ": Error, unknown averaging method..." << endl;
	// unknown method, fallback to simple
	SetAveragingMethod(SIMPLE, false);
}

void SAR_Calculation::SetNumLines(unsigned int numLines[3])
{
	Delete3DArray(m_Vx_Used,m_numLines);
	m_Vx_Used = NULL;
	Delete3DArray(m_Vx_Valid,m_numLines);
	m_Vx_Valid = NULL;

	for (int n=0;n<3;++n)
		m_numLines[n]=numLines[n];
}

void SAR_Calculation::SetCellWidth(float* cellWidth[3])
{
	for (int n=0;n<3;++n)
		m_cellWidth[n]=cellWidth[n];
}

float*** SAR_Calculation::CalcSAR(float*** SAR)
{
	if (CheckValid()==false)
	{
		cerr << "SAR_Calculation::CalcSAR: SAR calculation is invalid due to missing values... Abort..." << endl;
		return NULL;
	}
	if (m_avg_mass<=0)
		return CalcLocalSAR(SAR);
	return CalcAveragedSAR(SAR);
}

bool SAR_Calculation::CheckValid()
{
	for (int n=0;n<3;++n)
		if (m_cellWidth[n]==NULL)
			return false;
	if (m_E_field==NULL)
		return false;
	if ((m_J_field==NULL) && (m_cell_conductivity==NULL))
		return false;
	if (m_cell_density==NULL)
		return false;
	if (m_avg_mass<0)
		return false;
	return true;
}

double SAR_Calculation::CalcSARPower()
{
	if (CheckValid()==false)
	{
		cerr << "SAR_Calculation::CalcSARPower: SAR calculation is invalid due to missing values... Abort..." << endl;
		return 0;
	}
	double power=0;
	unsigned int pos[3];
	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
			{
				power += CalcLocalPowerDensity(pos)*CellVolume(pos);
			}
		}
	}
	return power;
}

double SAR_Calculation::CalcLocalPowerDensity(unsigned int pos[3])
{
	double l_pow=0;
	if (m_cell_conductivity==NULL)
	{
		l_pow  = abs(m_E_field[0][pos[0]][pos[1]][pos[2]]) * abs(m_J_field[0][pos[0]][pos[1]][pos[2]]);
		l_pow += abs(m_E_field[1][pos[0]][pos[1]][pos[2]]) * abs(m_J_field[1][pos[0]][pos[1]][pos[2]]);
		l_pow += abs(m_E_field[2][pos[0]][pos[1]][pos[2]]) * abs(m_J_field[2][pos[0]][pos[1]][pos[2]]);
	}
	else
	{
		l_pow  = m_cell_conductivity[pos[0]][pos[1]][pos[2]]*abs(m_E_field[0][pos[0]][pos[1]][pos[2]]) * abs(m_E_field[0][pos[0]][pos[1]][pos[2]]);
		l_pow += m_cell_conductivity[pos[0]][pos[1]][pos[2]]*abs(m_E_field[1][pos[0]][pos[1]][pos[2]]) * abs(m_E_field[1][pos[0]][pos[1]][pos[2]]);
		l_pow += m_cell_conductivity[pos[0]][pos[1]][pos[2]]*abs(m_E_field[2][pos[0]][pos[1]][pos[2]]) * abs(m_E_field[2][pos[0]][pos[1]][pos[2]]);
	}
	l_pow*=0.5;
	return l_pow;
}

float*** SAR_Calculation::CalcLocalSAR(float*** SAR)
{
	unsigned int pos[3];
	m_Valid=0;
	m_Used=0;
	m_Unused=0;
	m_AirVoxel=0;
	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
			{
				if (m_cell_density[pos[0]][pos[1]][pos[2]]>0)
				{
					++m_Valid;
					SAR[pos[0]][pos[1]][pos[2]]	= CalcLocalPowerDensity(pos)/m_cell_density[pos[0]][pos[1]][pos[2]];
				}
				else
				{
					++m_AirVoxel;
					SAR[pos[0]][pos[1]][pos[2]] = 0;
				}
			}
		}
	}
	return SAR;
}

int SAR_Calculation::FindFittingCubicalMass(unsigned int pos[3], float box_size, unsigned int start[3], unsigned int stop[3],
					float partial_start[3], float partial_stop[3], double &mass, double &volume, double &bg_ratio, int disabledFace, bool ignoreFaceValid)
{
	unsigned int mass_iterations = 0;
	double old_mass=0;
	double old_box_size=0;
	bool face_valid;
	bool mass_valid;
	bool voxel_valid;

	//iterate over cubical sizes to find fitting volume to mass
	while (mass_iterations<m_maxMassIterations)
	{
		// calculate cubical mass
		face_valid = GetCubicalMass(pos, box_size/2, start, stop, partial_start, partial_stop, mass, volume, bg_ratio, disabledFace);

		// check if found mass is valid
		mass_valid = abs(mass-m_avg_mass)<=m_massTolerance*m_avg_mass;
		voxel_valid = mass_valid && (face_valid==true) && (bg_ratio<m_maxBGRatio);

		if ((face_valid==false) && (mass<m_avg_mass*(1.0-m_massTolerance)) && (ignoreFaceValid==false))
		{
			// this is an invalid cube with a too small total mass --> increasing the box will not yield a valid cube
			return 1;
		}
		else if ((face_valid==false  || bg_ratio>=m_maxBGRatio) && (mass_valid))
		{
			// this is an invalid cube with a valid total mass
			if (ignoreFaceValid)
				return 0;
			return 2;
		}
		if (voxel_valid)
		{
			// valid cube found
			return 0;
		}

		// if no valid or finally invalid cube is found, calculate an alternaive cube size
		if (mass_iterations==0)
		{
			// on first interation, try a relative resize
			old_box_size=box_size;
			box_size*=pow(m_avg_mass/mass,1.0/3.0);
		}
		else
		{
			// on later interations, try a newton approach
			float new_box_size = box_size - (mass-m_avg_mass)/(mass-old_mass)*(box_size-old_box_size);
			old_box_size = box_size;
			box_size = new_box_size;
		}
		old_mass=mass;

		++mass_iterations;
	}

	// m_maxMassIterations iterations are exhausted...
	return -1;
}

bool SAR_Calculation::GetCubicalMass(unsigned int pos[3], double box_size, unsigned int start[3], unsigned int stop[3],
									 float partial_start[3], float partial_stop[3], double &mass, double &volume, double &bg_ratio, int disabledFace)
{
	if ((box_size<=0) || std::isnan(box_size) || std::isinf(box_size))
	{
		cerr << "SAR_Calculation::GetCubicalMass: critical error: invalid averaging box size!! EXIT" << endl;
		exit(-1);
	}
	bool face_valid=true;
	for (int n=0;n<3;++n)
	{
		// check start position
		start[n]=pos[n];
		partial_start[n]=1;
		float dist=m_cellWidth[n][pos[n]]/2;
		if (disabledFace==2*n)
			dist=box_size;
		else
		{
			while (dist<box_size)
			{
				// box is outside of domain
				if (start[n]==0)
				{
					partial_start[n]=-1;
					break;
				}
				--start[n];
				dist+=m_cellWidth[n][start[n]];
			}

			if (partial_start[n]!=-1)
			{   // box ends inside stop[n] voxel
				partial_start[n]=1-(dist-box_size)/m_cellWidth[n][start[n]];
			}
			if ((partial_start[n]!=-1) && (pos[n]==start[n]))
				partial_start[n]=2*box_size/m_cellWidth[n][start[n]];
		}

		// check stop position
		stop[n]=pos[n];
		partial_stop[n]=1;
		dist=m_cellWidth[n][pos[n]]/2;
		if (disabledFace==2*n+1)
			dist=box_size;
		else
		{
			while (dist<box_size)
			{
				// box is outside of domain
				if (stop[n]==m_numLines[n]-1)
				{
					partial_stop[n]=-1;
					break;
				}
				++stop[n];
				dist+=m_cellWidth[n][stop[n]];
			}

			if (partial_stop[n]!=-1)
			{   // box ends inside stop[n] voxel
				partial_stop[n]=1-(dist-box_size)/m_cellWidth[n][stop[n]];
			}
			if ((partial_stop[n]!=-1) && (pos[n]==stop[n]))
				partial_stop[n]=2*box_size/m_cellWidth[n][stop[n]];
		}
	}

	for (int n=0;n<3;++n)
	{
		if (partial_start[n]==-1)
			face_valid=false;
		if (partial_stop[n]==-1)
			face_valid=false;
	}

	mass = 0;
	volume = 0;
	double bg_volume=0;
	double weight[3];
	unsigned int f_pos[3];
	bool face_intersect[6] = {false,false,false,false,false,false};
	for (f_pos[0]=start[0];f_pos[0]<=stop[0];++f_pos[0])
	{
		weight[0]=1;
		if (f_pos[0]==start[0])
			weight[0]*=abs(partial_start[0]);
		if (f_pos[0]==stop[0])
			weight[0]*=abs(partial_stop[0]);

		for (f_pos[1]=start[1];f_pos[1]<=stop[1];++f_pos[1])
		{
			weight[1]=1;
			if (f_pos[1]==start[1])
				weight[1]*=abs(partial_start[1]);
			if (f_pos[1]==stop[1])
				weight[1]*=abs(partial_stop[1]);

			for (f_pos[2]=start[2];f_pos[2]<=stop[2];++f_pos[2])
			{
				weight[2]=1;
				if (f_pos[2]==start[2])
					weight[2]*=abs(partial_start[2]);
				if (f_pos[2]==stop[2])
					weight[2]*=abs(partial_stop[2]);

				mass += CellMass(f_pos)*weight[0]*weight[1]*weight[2];
				volume += CellVolume(f_pos)*weight[0]*weight[1]*weight[2];

				if (m_cell_density[f_pos[0]][f_pos[1]][f_pos[2]]==0)
					bg_volume += CellVolume(f_pos)*weight[0]*weight[1]*weight[2];
				else
				{
					for (int n=0;n<3;++n)
					{
						if (start[n]==f_pos[n])
							face_intersect[2*n]=true;
						if (stop[n]==f_pos[n])
							face_intersect[2*n+1]=true;
					}
				}
			}
		}
	}
	//check if all bounds have intersected a material boundary
	for (int n=0;n<6;++n)
		face_valid *= face_intersect[n];

	bg_ratio = bg_volume/volume;

	return face_valid;
}

float SAR_Calculation::CalcCubicalSAR(float*** SAR, unsigned int pos[3], unsigned int start[3], unsigned int stop[3], float partial_start[3], float partial_stop[3], bool assignUsed)
{
	double power_mass=0;
	double mass=0;
	double weight[3];
	unsigned int f_pos[3];
	for (f_pos[0]=start[0];f_pos[0]<=stop[0];++f_pos[0])
	{
		weight[0]=1;
		if (f_pos[0]==start[0])
			weight[0]*=abs(partial_start[0]);
		if (f_pos[0]==stop[0])
			weight[0]*=abs(partial_stop[0]);

		for (f_pos[1]=start[1];f_pos[1]<=stop[1];++f_pos[1])
		{
			weight[1]=1;
			if (f_pos[1]==start[1])
				weight[1]*=abs(partial_start[1]);
			if (f_pos[1]==stop[1])
				weight[1]*=abs(partial_stop[1]);

			for (f_pos[2]=start[2];f_pos[2]<=stop[2];++f_pos[2])
			{
				weight[2]=1;
				if (f_pos[2]==start[2])
					weight[2]*=abs(partial_start[2]);
				if (f_pos[2]==stop[2])
					weight[2]*=abs(partial_stop[2]);

				if (m_cell_density[f_pos[0]][f_pos[1]][f_pos[2]]>=0)
				{
					mass += CellMass(f_pos)*weight[0]*weight[1]*weight[2];
					power_mass += CalcLocalPowerDensity(f_pos) * CellVolume(f_pos)*weight[0]*weight[1]*weight[2];
				}
			}
		}
	}
	float vx_SAR = power_mass/mass;
	if (SAR==NULL)
		return vx_SAR;

	SAR[pos[0]][pos[1]][pos[2]]=vx_SAR;

	if (assignUsed==false)
		return vx_SAR;

	// assign SAR to all used voxel
	bool is_partial[3];
	for (f_pos[0]=start[0];f_pos[0]<=stop[0];++f_pos[0])
	{
		if ( ((f_pos[0]==start[0]) && (partial_start[0]!=1)) || ((f_pos[0]==stop[0]) && (partial_stop[0]!=1)) )
			is_partial[0]=true;
		else
			is_partial[0]=false;

		for (f_pos[1]=start[1];f_pos[1]<=stop[1];++f_pos[1])
		{
			if ( ((f_pos[1]==start[1]) && (partial_start[1]!=1)) || ((f_pos[1]==stop[1]) && (partial_stop[1]!=1)) )
				is_partial[1]=true;
			else
				is_partial[1]=false;

			for (f_pos[2]=start[2];f_pos[2]<=stop[2];++f_pos[2])
			{
				if ( ((f_pos[2]==start[2]) && (partial_start[2]!=1)) || ((f_pos[2]==stop[2]) && (partial_stop[2]!=1)) )
					is_partial[2]=true;
				else
					is_partial[2]=false;

				if ( (!is_partial[0] && !is_partial[1] && !is_partial[2]) || m_markPartialAsUsed)
				{
					if (!m_Vx_Valid[f_pos[0]][f_pos[1]][f_pos[2]] && (m_cell_density[f_pos[0]][f_pos[1]][f_pos[2]]>0))
					{
						m_Vx_Used[f_pos[0]][f_pos[1]][f_pos[2]]=true;
						SAR[f_pos[0]][f_pos[1]][f_pos[2]]=max(SAR[f_pos[0]][f_pos[1]][f_pos[2]], vx_SAR);
					}
				}
			}
		}
	}
	return vx_SAR;
}


float*** SAR_Calculation::CalcAveragedSAR(float*** SAR)
{
	unsigned int pos[3];
	m_Vx_Used = Create3DArray<bool>(m_numLines);
	m_Vx_Valid = Create3DArray<bool>(m_numLines);

	double voxel_volume;
	double total_mass;
	unsigned int start[3];
	unsigned int stop[3];
	float partial_start[3];
	float partial_stop[3];
	double bg_ratio;
	int EC=0;

	// debug counter
	unsigned int cnt_case1=0;
	unsigned int cnt_case2=0;
	unsigned int cnt_NoConvergence=0;

	m_Valid=0;
	m_Used=0;
	m_Unused=0;
	m_AirVoxel=0;

	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
			{
				if (m_cell_density[pos[0]][pos[1]][pos[2]]==0)
				{
					SAR[pos[0]][pos[1]][pos[2]] = 0;
					++m_AirVoxel;
					continue;
				}

				// guess an initial box size and find a fitting cube
				EC = FindFittingCubicalMass(pos, pow(m_avg_mass/m_cell_density[pos[0]][pos[1]][pos[2]],1.0/3.0)/2, start, stop,
											partial_start, partial_stop, total_mass, voxel_volume, bg_ratio, -1, m_IgnoreFaceValid);

				if (EC==0)
				{
					m_Vx_Valid[pos[0]][pos[1]][pos[2]] = true;
					m_Vx_Used[pos[0]][pos[1]][pos[2]] = true;
					++m_Valid;
					CalcCubicalSAR(SAR, pos, start, stop, partial_start, partial_stop, true);
				}
				else if (EC==1)
					++cnt_case1;
				else if (EC==2)
					++cnt_case2;
				else if (EC==-1)
					++cnt_NoConvergence;
				else
					cerr << "other EC" << EC << endl;
			}
		}
	}
	if (cnt_NoConvergence>0)
	{
		cerr << "SAR_Calculation::CalcAveragedSAR: Warning, for some voxel a valid averaging cube could not be found (no convergence)... " << endl;
	}
	if (m_DebugLevel>0)
	{
		cerr << "Number of invalid cubes (case 1): " << cnt_case1 << endl;
		cerr << "Number of invalid cubes (case 2): " << cnt_case2 << endl;
		cerr << "Number of invalid cubes (failed to converge): " << cnt_NoConvergence << endl;
	}

	// count all used and unused etc. + special handling of unused voxels!!
	m_Used=0;
	m_Unused=0;
	for (pos[0]=0;pos[0]<m_numLines[0];++pos[0])
	{
		for (pos[1]=0;pos[1]<m_numLines[1];++pos[1])
		{
			for (pos[2]=0;pos[2]<m_numLines[2];++pos[2])
			{
				if (!m_Vx_Valid[pos[0]][pos[1]][pos[2]] && m_Vx_Used[pos[0]][pos[1]][pos[2]])
					++m_Used;
				if ((m_cell_density[pos[0]][pos[1]][pos[2]]>0) && !m_Vx_Valid[pos[0]][pos[1]][pos[2]] && !m_Vx_Used[pos[0]][pos[1]][pos[2]])
				{
					++m_Unused;

					SAR[pos[0]][pos[1]][pos[2]] = 0;
					double unused_volumes[6];
					float unused_SAR[6];

					double min_Vol=FLT_MAX;

					// special handling of unused voxels:
					for (int n=0;n<6;++n)
					{
						EC = FindFittingCubicalMass(pos, pow(m_avg_mass/m_cell_density[pos[0]][pos[1]][pos[2]],1.0/3.0)/2, start, stop,
													partial_start, partial_stop, total_mass, unused_volumes[n], bg_ratio, n, true);
						if ((EC!=0) && (EC!=2))
						{
							// this should not happen
							cerr << "SAR_Calculation::CalcAveragedSAR: Error handling unused voxels, can't find fitting cubical averaging volume' " << endl;
							unused_SAR[n]=0;
						}
						else
						{
							unused_SAR[n]=CalcCubicalSAR(NULL, pos, start, stop, partial_start, partial_stop, false);
							min_Vol = min(min_Vol,unused_volumes[n]);
						}
					}
					for (int n=0;n<6;++n)
					{
						if (unused_volumes[n]<=m_UnusedRelativeVolLimit*min_Vol)
							SAR[pos[0]][pos[1]][pos[2]] = max(SAR[pos[0]][pos[1]][pos[2]],unused_SAR[n]);
					}
				}
			}
		}
	}

	if (m_Valid+m_Used+m_Unused+m_AirVoxel!=m_numLines[0]*m_numLines[1]*m_numLines[2])
	{
		cerr << "SAR_Calculation::CalcAveragedSAR: critical error, mismatch in voxel status count... EXIT" << endl;
		exit(1);
	}

	if (m_DebugLevel>0)
		cerr << "SAR_Calculation::CalcAveragedSAR: Stats: Valid=" << m_Valid << " Used=" << m_Used << " Unused=" << m_Unused << " Air-Voxel=" << m_AirVoxel << endl;

	return SAR;
}

double SAR_Calculation::CellVolume(unsigned int pos[3])
{
	if (m_cell_volume)
		return m_cell_volume[pos[0]][pos[1]][pos[2]];

	double volume=1;
	for (int n=0;n<3;++n)
		volume*=m_cellWidth[n][pos[n]];
	return volume;
}

double SAR_Calculation::CellMass(unsigned int pos[3])
{
	return m_cell_density[pos[0]][pos[1]][pos[2]]*CellVolume(pos);
}

