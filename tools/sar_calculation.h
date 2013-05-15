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

#ifndef SAR_CALCULATION_H
#define SAR_CALCULATION_H

#include <complex>

class SAR_Calculation
{
public:
	SAR_Calculation();

	enum SARAveragingMethod { IEEE_C95_3, IEEE_62704, SIMPLE};

	//! Reset and initialize all values (will keep all SAR settings)
	void Reset();

	//! Set the debug level
	void SetDebugLevel(int level) {m_DebugLevel=level;}

	//! Set the used averaging method
	void SetAveragingMethod(SARAveragingMethod method, bool silent=false);

	//! Set the used averaging method
	void SetAveragingMethod(std::string method, bool silent=false);

	//! Set number of lines in all direcitions. (mandatory information)
	void SetNumLines(unsigned int numLines[3]);
	//! Set cell width in all direcitions. (mandatory information for averaging)
	void SetCellWidth(float* cellWidth[3]);

	//! Set the averaging mash. (mandatory information for averaging)
	void SetAveragingMass(float mass) {m_avg_mass=mass;}

	//! Set the cell volumes (optional for speedup)
	void SetCellVolumes(float*** cell_volume) {m_cell_volume=cell_volume;}

	//! Set the cell densities (mandatory information)
	void SetCellDensities(float*** cell_density) {m_cell_density=cell_density;}

	//! Set the cell conductivities (mandatory if no current density field is given)
	void SetCellCondictivity(float*** cell_conductivity) {m_cell_conductivity=cell_conductivity;}

	//! Set the electric field (mandatory information)
	void SetEField(std::complex<float>**** field) {m_E_field=field;}
	//! Set the current density field (mandatory if no conductivity distribution is given)
	void SetJField(std::complex<float>**** field) {m_J_field=field;}

	//! Calculate the SAR, requires a preallocated 3D array
	float*** CalcSAR(float*** SAR);

	//! Calculate the total power dumped
	double CalcSARPower();

protected:
	unsigned int m_numLines[3];
	float* m_cellWidth[3];

	float m_avg_mass;
	float*** m_cell_volume;
	float*** m_cell_density;
	float*** m_cell_conductivity;
	std::complex<float>**** m_E_field;
	std::complex<float>**** m_J_field;

	bool*** m_Vx_Used;
	bool*** m_Vx_Valid;

	unsigned int m_Valid;
	unsigned int m_Used;
	unsigned int m_Unused;
	unsigned int m_AirVoxel;

	int m_DebugLevel;

	/*********** SAR calculation parameter and settings ***********/
	float m_massTolerance;
	unsigned int m_maxMassIterations;
	float m_maxBGRatio;
	bool m_markPartialAsUsed;
	float m_UnusedRelativeVolLimit;
	bool m_IgnoreFaceValid;

	/*********** SAR calculations methods ********/
	double CalcLocalPowerDensity(unsigned int pos[3]);

	//! Calculate the local SAR
	float*** CalcLocalSAR(float*** SAR);

	/****** start SAR averaging and all necessary methods ********/
	//! Calculate the averaged SAR
	float*** CalcAveragedSAR(float*** SAR);

	int FindFittingCubicalMass(unsigned int pos[3], float box_size, unsigned int start[3], unsigned int stop[3],
						float partial_start[3], float partial_stop[3], double &mass, double &volume, double &bg_ratio, int disabledFace=-1, bool ignoreFaceValid=false);
	bool GetCubicalMass(unsigned int pos[3], double box_size, unsigned int start[3], unsigned int stop[3],
						float partial_start[3], float partial_stop[3], double &mass, double &volume, double &bg_ratio, int disabledFace=-1);

	float CalcCubicalSAR(float*** SAR, unsigned int pos[3], unsigned int start[3], unsigned int stop[3], float partial_start[3], float partial_stop[3], bool assignUsed=false);
	/****** end SAR averaging and all necessary methods ********/

	bool CheckValid();
	double CellVolume(unsigned int pos[3]);
	double CellMass(unsigned int pos[3]);
};

#endif // SAR_CALCULATION_H
