/*
*	Copyright (C) 2012-2025 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include <vector>
#include <mutex>
#include <atomic>
#include "arraylib/array_ijk.h"
#include "arraylib/array_nijk.h"

#if defined(WIN32)
	#ifdef BUILD_OPENEMS_LIB
	#define SAR_EXPORT __declspec(dllexport)
	#else
	#define SAR_EXPORT __declspec(dllimport)
	#endif
#else
#define SAR_EXPORT
#endif

class HDF5_File_Writer;

class SAR_EXPORT SAR_Calculation
{
public:
	SAR_Calculation();
	~SAR_Calculation();

	enum SARAveragingMethod { IEEE_C95_3, IEEE_62704, SIMPLE};

	//! Reset and initialize all values (will keep all SAR settings)
	void Reset();

	//! Set the debug level
	void SetDebugLevel(int level) {m_DebugLevel=level;}

	//! Enable progress output
	void EnableProgress(bool enable=true) {m_progress=enable;}

	//! Set the used averaging method
	bool SetAveragingMethod(SARAveragingMethod method, bool silent=false);

	//! Set the used averaging method
	bool SetAveragingMethod(std::string method, bool silent=false);

	//! Set number of lines in all direcitions. (mandatory information)
	void SetNumLines(unsigned int numLines[3]);

	//! Define a sub range for SAR calculation for the given direction. Number of lines must already be set!
	void SetSubRange(int ny, unsigned int start, unsigned int stop);

	//! Set cell width in all direcitions. (mandatory information for averaging)
	void SetCellWidth(double* cellWidth[3]);

	//! Set the averaging mash. (mandatory information for averaging)
	void SetAveragingMass(float mass) {m_avg_mass=mass;}

	//! Set the cell volumes (optional for speedup)
	void SetCellVolumes(ArrayLib::ArrayIJK<float>* cell_volume) {m_cell_volume=cell_volume;}

	//! Set the cell densities (mandatory information)
	void SetCellDensities(ArrayLib::ArrayIJK<float>* cell_density) {m_cell_density=cell_density;}

	//! Add an electric field
	void AddEFieldAndCondictivity(float freq, ArrayLib::ArrayNIJK<std::complex<float>>* e_field, ArrayLib::ArrayIJK<float>* cell_conductivity);

	//! Add an electric field and current density field
	void AddEFieldAndJField(float freq, ArrayLib::ArrayNIJK<std::complex<float>>* e_field, ArrayLib::ArrayNIJK<std::complex<float>>* j_field);

	//! Enable the recording of special averaging cube stats
	void EnableCubeStats() {m_record_cube_stats=true;}

	//! Limit the SAR calculation area to a box with the highest local power densities
	void EnableAutoRange(double dBmax) {m_autoRange=dBmax;}

	//! Calculate the SAR, requires a preallocated 3D array
	bool CalcSAR(unsigned int numThreads=0);

	//! Get the total power dumped for frequency index n
	double GetSARPower(size_t freq_n) {return m_power.at(freq_n);}

	//! Get the SAR results for frequency index n
	ArrayLib::ArrayIJK<float>* GetSAR(size_t freq_n) {return m_SAR.at(freq_n);}

	//! Calculate the total mass
	double CalcTotalMass();

	//! Read all raw data from an hdf5 file
	bool ReadFromHDF5(std::string h5_fn);

	//! Create and write an hdf5 output file (incl. mesh etc.)
	bool WriteToHDF5(std::string out_name, bool legacyHDF5=false);

	//! Write the results to an hdf5 output file (mesh information and some other details are not included here)
	bool WriteToHDF5(HDF5_File_Writer &out_file, bool legacyHDF5=false);

	//! Read and write from an hdf5 file
	bool CalcFromHDF5(std::string h5_fn, std::string out_name, bool legacyHDF5=false, unsigned int numThreads=0);

protected:
	// raw input data
	unsigned int m_numLines[3];
	double* m_lines[3]; // may not be set in all cases

	bool m_cleanup_cell_data = false;
	double* m_cellWidth[3];
	ArrayLib::ArrayIJK<float>* m_cell_volume;
	ArrayLib::ArrayIJK<float>* m_cell_density;

	std::vector<float> m_freq;
	std::vector<ArrayLib::ArrayIJK<float>*> m_local_cell_power_density; // precalculated local cell power density
	std::vector<float> m_local_cell_max_power_density; // precalculated max local cell power density (e.g. needed for auto range)

	// save some statistical data
	size_t m_Valid;
	size_t m_Used;
	size_t m_Unused;
	size_t m_AirVoxel;

	// output data
	ArrayLib::ArrayIJK<unsigned int> m_posInOuput;
	std::vector<unsigned int> m_cellIndicies[3];  // list of output indecies into the raw data mesh
	ArrayLib::ArrayIJK<unsigned char> m_cube_type;
	ArrayLib::ArrayIJK<float> m_cube_mass;
	ArrayLib::ArrayIJK<float> m_cube_volume;

	std::vector<double> m_power;
	std::vector<float> m_maxSAR;
	std::vector<std::array<unsigned int, 3>> m_maxSAR_Idx;
	std::vector<ArrayLib::ArrayIJK<float>*> m_SAR;

	int m_DebugLevel;
	double m_duration;
	double m_autoRange = 0;
	bool m_record_cube_stats = false;
	bool m_progress = false;

	/*********** SAR calculation parameter and settings ***********/
	float m_avg_mass;
	float m_massTolerance;
	unsigned int m_maxMassIterations;
	float m_maxBGRatio;
	bool m_markPartialAsUsed;
	float m_UnusedRelativeVolLimit;
	bool m_IgnoreFaceValid;

	/*********** SAR calculations methods ********/
	void DoAutoRange();
	void InitSAR();

	//! Calculate the local SAR
	bool CalcLocalSAR();

	/****** start SAR averaging and all necessary methods ********/
	//! Calculate the averaged SAR
	unsigned int m_numThreads=0;
	std::atomic<size_t> m_progressCounter{0};
	std::mutex m_resultMutex;
	std::mutex m_planeLocks[128];
	std::atomic<size_t> m_nextStepChunk{0};
	size_t m_step1_case1, m_step1_case2, m_step1_no_conv;
	bool CalcAvgStep1SAR(unsigned int t_id, ArrayLib::ArrayIJK<bool> &Vx_Valid, ArrayLib::ArrayIJK<bool> &Vx_Used);
	bool CalcAvgStep2SAR(unsigned int t_id, ArrayLib::ArrayIJK<bool> &Vx_Valid, ArrayLib::ArrayIJK<bool> &Vx_Used);
	bool CalcAveragedSAR(unsigned int numThreads=0);

	int FindFittingCubicalMass(unsigned int pos[3], float box_size, unsigned int start[3], unsigned int stop[3],
						float partial_start[3], float partial_stop[3], double &mass, double &volume, double &bg_ratio, int disabledFace=-1, bool ignoreFaceValid=false);
	bool GetCubicalMass(unsigned int pos[3], double box_size, unsigned int start[3], unsigned int stop[3],
						float partial_start[3], float partial_stop[3], double &mass, double &volume, double &bg_ratio, int disabledFace=-1);

	void CalcCubicalSAR(double* vx_sar, unsigned int start[3], unsigned int stop[3],
						 float partial_start[3], float partial_stop[3]);
	void AssignUsedSAR(double* vx_sar, unsigned int start[3], unsigned int stop[3],
						 float partial_start[3], float partial_stop[3],
						 ArrayLib::ArrayIJK<bool> &Vx_Used, ArrayLib::ArrayIJK<bool> &Vx_Valid);
	/****** end SAR averaging and all necessary methods ********/

	bool CheckValid();
};

#endif // SAR_CALCULATION_H
