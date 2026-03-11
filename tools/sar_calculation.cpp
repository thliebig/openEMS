/*
*	Copyright (C) 2012-2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include <chrono>
#include "sar_calculation.h"
#include "cfloat"
#include "useful.h"
#include "array_ops.h"
#include "global.h"
#include "hdf5_file_reader.h"

#include "arraylib/array_nijk.h"

using namespace std;

SAR_Calculation::SAR_Calculation()
{
	m_DebugLevel = 0;
	m_avg_mass = 0;
	SetAveragingMethod(SIMPLE, true);
	m_cellWidth[0]=m_cellWidth[1]=m_cellWidth[2]=NULL;
	m_lines[0]=m_lines[1]=m_lines[2]=NULL;
	Reset();
}

SAR_Calculation::~SAR_Calculation()
{
	Reset();
}

void SAR_Calculation::Reset()
{
	m_duration = 0;
	for (int n=0;n<3;++n)
	{
		m_cellIndicies[n].clear();
		m_numLines[n]==0;
		if (m_cleanup_cell_data)
		{
			delete[] m_cellWidth[n];
			m_cellWidth[n]=NULL;
		}
		delete[] m_lines[n];
		m_lines[n]=NULL;
	}
	m_posInOuput.Reset();
	if (m_cleanup_cell_data)
	{
		delete m_cell_volume;
		delete m_cell_density;
	}
	m_cell_volume = NULL;
	m_cell_density = NULL;
	for (size_t i=0;i<m_local_cell_power_density.size();++i)
		delete m_local_cell_power_density.at(i);
	m_local_cell_power_density.clear();

	for (size_t i=0;i<m_SAR.size();++i)
		delete m_SAR.at(i);
	m_SAR.clear();
}

bool SAR_Calculation::SetAveragingMethod(string method, bool silent)
{
	if (method.compare("IEEE_C95_3")==0)
		return SetAveragingMethod(IEEE_C95_3, silent);
	if (method.compare("IEEE_62704")==0)
		return SetAveragingMethod(IEEE_62704, silent);
	if ((method.compare("Simple")==0) || (method.compare("SIMPLE")==0))
		return SetAveragingMethod(SIMPLE, silent);

	cerr << __func__ << ": Error, """ << method << """ is an unknown averaging method..." << endl;
	return false;
}

bool SAR_Calculation::SetAveragingMethod(SARAveragingMethod method, bool silent)
{
	if (method==IEEE_62704)
	{
		m_massTolerance = 0.000001;  // yes, this low mass tolerance is in IEC/IEEE 62704-1
		m_maxMassIterations = 100;
		m_maxBGRatio = 0.1;
		m_markPartialAsUsed = false;
		m_UnusedRelativeVolLimit = 1.05;
		m_IgnoreFaceValid = false;
		if (!silent)
			cerr << __func__ << ": Setting averaging method to IEEE_62704" << endl;
		return true;
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
		return true;
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
		return true;
	}

	cerr << __func__ << ": Error, unknown averaging method..." << endl;
	return false;
}

void SAR_Calculation::SetNumLines(unsigned int numLines[3])
{
	for (int n=0;n<3;++n)
	{
		m_numLines[n]=numLines[n];
		m_cellIndicies[n].clear();
		m_cellIndicies[n].reserve(m_numLines[n]);
		for (unsigned int i=0;i<m_numLines[n];++i)
			m_cellIndicies[n].push_back(i);
	}
}

void SAR_Calculation::SetSubRange(int ny, unsigned int start, unsigned int stop)
{
	if ((start>=stop) || (start>=m_numLines[ny]))
	{
		cerr << __func__ << " Invalid range!" << endl;
		return;
	}

	stop = min(stop, m_numLines[ny]);
	m_cellIndicies[ny].clear();
	m_cellIndicies[ny].reserve(stop-start);
	if (m_DebugLevel>0)
		cout << "SetSubRange: Dir " << ny << " Start " << start << " Stop " << stop << endl;
	for (unsigned int i=start;i<stop;++i)
		m_cellIndicies[ny].push_back(i);
}

void SAR_Calculation::SetCellWidth(double* cellWidth[3])
{
	for (int n=0;n<3;++n)
		m_cellWidth[n]=cellWidth[n];
}

void SAR_Calculation::AddEFieldAndCondictivity(float freq, ArrayLib::ArrayNIJK<std::complex<float>>* e_field, ArrayLib::ArrayIJK<float>* cell_conductivity)
{
	ArrayLib::ArrayIJK<float>* loc_cpd = new ArrayLib::ArrayIJK<float>("loc_cpd", m_numLines);
	float* loc_cpd_data = loc_cpd->data();
	float* cell_vol = m_cell_volume->data();
	float* cell_dens = m_cell_density->data();
	double max_lp=0;
	unsigned int loc=0;
	unsigned int pos[3];
	double l_pow=0;
	double power=0;
	if (e_field->extent(0)!=3)
		throw;
	for (int n=0;n<3;++n)
		if (e_field->extent(n+1) != m_numLines[n])
			throw;

	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
			{
				l_pow  = (*cell_conductivity)(pos[0], pos[1], pos[2])*abs((*e_field)(0, pos[0], pos[1], pos[2])) * abs((*e_field)(0, pos[0], pos[1], pos[2]));
				l_pow += (*cell_conductivity)(pos[0], pos[1], pos[2])*abs((*e_field)(1, pos[0], pos[1], pos[2])) * abs((*e_field)(1, pos[0], pos[1], pos[2]));
				l_pow += (*cell_conductivity)(pos[0], pos[1], pos[2])*abs((*e_field)(2, pos[0], pos[1], pos[2])) * abs((*e_field)(2, pos[0], pos[1], pos[2]));
				l_pow *= 0.5*(cell_dens[loc]>0);
				max_lp = max(max_lp, l_pow);
				power += l_pow*cell_vol[loc];
				loc_cpd_data[loc] = l_pow;
				++loc;
			}
		}
	}

	m_local_cell_power_density.push_back(loc_cpd);
	m_freq.push_back(freq);
	m_power.push_back(power);
	m_local_cell_max_power_density.push_back(max_lp);
}

void SAR_Calculation::AddEFieldAndJField(float freq, ArrayLib::ArrayNIJK<std::complex<float>>* e_field, ArrayLib::ArrayNIJK<std::complex<float>>* j_field)
{
	ArrayLib::ArrayIJK<float>* loc_cpd = new ArrayLib::ArrayIJK<float>("loc_cpd", m_numLines);
	float* loc_cpd_data = loc_cpd->data();
	float* cell_vol = m_cell_volume->data();
	float* cell_dens = m_cell_density->data();
	double max_lp=0;
	unsigned int loc;
	unsigned int pos[3];
	double l_pow;
	double power=0;
	for (pos[0]=0; pos[0]<m_numLines[0]; ++pos[0])
	{
		for (pos[1]=0; pos[1]<m_numLines[1]; ++pos[1])
		{
			for (pos[2]=0; pos[2]<m_numLines[2]; ++pos[2])
			{
				l_pow  = abs((*e_field)(0, pos[0], pos[1], pos[2])) * abs((*j_field)(0, pos[0], pos[1], pos[2]));
				l_pow += abs((*e_field)(1, pos[0], pos[1], pos[2])) * abs((*j_field)(1, pos[0], pos[1], pos[2]));
				l_pow += abs((*e_field)(2, pos[0], pos[1], pos[2])) * abs((*j_field)(2, pos[0], pos[1], pos[2]));
				l_pow *= 0.5*(cell_dens[loc]>0);
				max_lp = max(max_lp, l_pow);
				power += l_pow*cell_vol[loc];
				loc_cpd_data[loc] = l_pow;
				++loc;
			}
		}
	}
	m_local_cell_power_density.push_back(loc_cpd);
	m_freq.push_back(freq);
	m_power.push_back(power);
	m_local_cell_max_power_density.push_back(max_lp);
}

bool SAR_Calculation::CalcSAR()
{
	if (CheckValid()==false)
	{
		cerr << "SAR_Calculation::CalcSAR: SAR calculation is invalid due to missing values... Abort..." << endl;
		return false;
	}

	auto start = std::chrono::high_resolution_clock::now();

	InitSAR();

	bool rc = false;
	if (m_avg_mass<=0)
		rc = CalcLocalSAR();
	else
		rc = CalcAveragedSAR();

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> duration = end - start;
	m_duration = duration.count()/1000;
	if (m_DebugLevel>0)
		cout << "Processing Time: " << m_duration << "s" << endl;

	return rc;
}

bool SAR_Calculation::CheckValid()
{
	for (int n=0;n<3;++n)
		if (m_cellWidth[n]==NULL)
			return false;

	if (m_local_cell_power_density.size()==0)
		return false;
	if ((m_cell_density==NULL) || (!m_cell_density->valid()))
		return false;
	if ((m_cell_volume==NULL) || (!m_cell_volume->valid()))
		return false;
	if (m_avg_mass<0)
		return false;

	// vaildate sizes
	for (int n=0;n<3;++n)
	{
		if (m_numLines[n]==0)
			return false;
		if (m_cell_density->extent(n) != m_numLines[n])
			return false;
		if (m_cell_volume->extent(n) != m_numLines[n])
			return false;
	}
	return true;
}

double SAR_Calculation::CalcTotalMass()
{
	double mass = 0;
	float* cell_density = m_cell_density->data();
	float* cell_volume = m_cell_volume->data();
	for (unsigned int n=0;n<m_cell_density->size();++n)
		mass += cell_density[n]*cell_volume[n];
	return mass;
}

bool SAR_Calculation::ReadFromHDF5(std::string h5_fn)
{
	Reset();
	if (m_DebugLevel>0)
	{
		cout << "Input File: " << h5_fn << endl;
		cout << "Averaging mass: " << m_avg_mass << endl;
	}
	HDF5_File_Reader h5_reader(h5_fn);
	if (!h5_reader.IsValid())
	{
		cerr << "SAR_Calculation::CalcFromHDF5: HDF5 file is invalid... Abort..." << endl;
		Reset();
		return false;
	}

	vector<float> freq;
	if ((h5_reader.ReadFrequencies(freq)==false) || (freq.size()==0))
	{
		cerr << "SAR_Calculation::CalcFromHDF5: No frequencies found..." << endl;
		Reset();
		return false;
	}
	if (m_DebugLevel>0)
		cout << "SAR_Calculation::CalcFromHDF5: Total number of frequencies: " << freq.size() << endl;

	// read mesh
	// double* lines[3]={NULL,NULL,NULL};
	// unsigned int numLines[3];
	int meshType;
	if (h5_reader.ReadMesh(m_lines, m_numLines, meshType) == false)
	{
		cerr << "SAR_Calculation::CalcFromHDF5: Error reading  mesh..." << endl;
		Reset();
		return false;
	}
	if (m_DebugLevel>0)
		cout << "SAR_Calculation::CalcFromHDF5: Mesh read complete. Total number of cells: " << m_numLines[0]*m_numLines[1]*m_numLines[2] << endl;

	// read cell width
	double* cellWidth[3]={NULL,NULL,NULL};
	unsigned int cw_numLines[3];
	int cw_meshType;
	if (h5_reader.ReadMesh(cellWidth, cw_numLines, cw_meshType, "/CellWidth") == false)
	{
		cerr << "SAR_Calculation::CalcFromHDF5: Error reading cell width..." << endl;
		Reset();
		return false;
	}

	m_cleanup_cell_data = true;
	m_cell_density = new ArrayLib::ArrayIJK<float>();
	if (!h5_reader.ReadScalarDataSet<float>("/CellData/Density", *m_cell_density))
	{
		cerr << "SAR_Calculation::CalcFromHDF5: Error reading  Density..." << endl;
		Reset();
		return false;
	}

	ArrayLib::ArrayIJK<float> cell_conductivity;
	if (!h5_reader.ReadScalarDataSet<float>("/CellData/Conductivity", cell_conductivity))
	{
		cerr << "SAR_Calculation::CalcFromHDF5: Error reading conductivity..." << endl;
		Reset();
		return false;
	}

	m_cell_volume = new ArrayLib::ArrayIJK<float>();
	if (!h5_reader.ReadScalarDataSet<float>("/CellData/Volume", *m_cell_volume))
	{
		cerr << "SAR_Calculation::CalcFromHDF5: Error reading cell volumes..." << endl;
		Reset();
		return false;
	}

	SetNumLines(m_numLines);
	SetCellWidth(cellWidth);

	ArrayLib::ArrayNIJK<complex<float>> E_fd_data;
	for (size_t n=0;n<freq.size();++n)
	{
		if (m_DebugLevel>0)
			cout << "SAR_Calculation::CalcFromHDF5: Read raw data for freq-index: " << n << endl;
		//E_fd_data = h5_reader.GetFDVectorData(n,data_size);
		if (!h5_reader.GetFDVectorData(n, E_fd_data) || (E_fd_data.extent(1)!=m_numLines[0]) || (E_fd_data.extent(2)!=m_numLines[1]) || (E_fd_data.extent(3)!=m_numLines[2]) )
		{
			cerr << E_fd_data.extent(1) << "," << E_fd_data.extent(2) << "," << E_fd_data.extent(3) << endl;
			cerr << "SAR_Calculation::CalcFromHDF5: FD data size mismatch... " << endl;
			Reset();
			return false;
		}
		AddEFieldAndCondictivity(freq.at(n), &E_fd_data, &cell_conductivity);
	}
	return true;
}

bool SAR_Calculation::WriteToHDF5(std::string out_name, bool legacyHDF5)
{
	if (m_DebugLevel>0)
		cout << "Output File: " << out_name << endl;

	HDF5_File_Writer out_file(out_name);

	unsigned int out_num_lines[3] = {(unsigned int)m_cellIndicies[0].size(),(unsigned int)m_cellIndicies[1].size(),(unsigned int)m_cellIndicies[2].size()};
	double* out_lines[3]={NULL,NULL,NULL};
	for (int n=0;n<3;++n)
	{
		out_lines[n] = new double[out_num_lines[n]];
		for (unsigned int i=0;i<out_num_lines[n];++i)
			out_lines[n][i] = m_lines[n][m_cellIndicies[n].at(i)];
	}

	out_file.WriteRectMesh(out_num_lines, out_lines, 0, 1);
	for (int n=0;n<3;++n)
		delete[] out_lines[n];

	out_file.WriteAtrribute("/","openEMS_HDF5_version", OPENEMS_HDF5_VERSION);
	out_file.SetCurrentGroup("/FieldData/FD");
	out_file.WriteAtrribute("/FieldData/FD","frequency",m_freq);

	return WriteToHDF5(out_file, legacyHDF5);
}

bool SAR_Calculation::WriteToHDF5(HDF5_File_Writer &out_file, bool legacyHDF5)
{
	size_t dims[3];// = new size_t[3];
	for (size_t n=0;n<3;++n)
		dims[n]=m_cellIndicies[n].size();

	double mass = CalcTotalMass();
	if (m_DebugLevel>0)
		cout << "Total Mass: " << mass << endl;
	out_file.WriteAtrribute("/","mass",mass);

	for (size_t n=0;n<m_freq.size();++n)
	{

		// CalcSAR(SAR);
		float pwr = GetSARPower(n);
		if (m_DebugLevel>0)
		{
			cout << "********" << endl;
			cout << "Frequency: " << m_freq.at(n) << endl;
			cout << "Total Power: " << pwr << endl;
			cout << "Avg SAR: " << pwr/mass << endl;
			cout << "max. SAR: " << m_maxSAR.at(n) << endl;
			cout << "max. SAR Pos.Idx.: " << m_maxSAR_Idx.at(n)[0] << ", " << m_maxSAR_Idx.at(n)[1] << ", " << m_maxSAR_Idx.at(n)[2]  << endl;
		}

		stringstream ss;
		ss << "f" << n;
		if (out_file.WriteScalarField<float>(ss.str(), *m_SAR.at(n), legacyHDF5)==false)
			cerr << "SAR_Calculation::CalcFromHDF5: can't dump to file...! " << endl;
		if (m_cube_type.data()!=NULL)
			if (out_file.WriteData(ss.str() + "_CubeType", H5T_NATIVE_UCHAR, m_cube_type.data(), 3, dims)==false)
				cerr << "SAR_Calculation::CalcFromHDF5: can't dump to file...! " << endl;
		if (m_cube_mass.data()!=NULL)
			if (out_file.WriteData(ss.str() + "_CubeMass", H5T_NATIVE_FLOAT, m_cube_mass.data(), 3, dims)==false)
				cerr << "SAR_Calculation::CalcFromHDF5: can't dump to file...! " << endl;
		if (m_cube_volume.data()!=NULL)
			if (out_file.WriteData(ss.str() + "_CubeVol", H5T_NATIVE_FLOAT, m_cube_volume.data(), 3, dims)==false)
				cerr << "SAR_Calculation::CalcFromHDF5: can't dump to file...! " << endl;
		if (out_file.WriteAtrribute("/FieldData/FD/"+ss.str(),"frequency",m_freq.at(n))==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
		if (out_file.WriteAtrribute("/FieldData/FD/"+ss.str(),"power",pwr)==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
		if (out_file.WriteAtrribute("/FieldData/FD/"+ss.str(),"maxSAR", m_maxSAR.at(n))==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
		std::vector<unsigned int> idx = {m_maxSAR_Idx.at(n)[0],m_maxSAR_Idx.at(n)[1],m_maxSAR_Idx.at(n)[2]};
		if (out_file.WriteAtrribute("/FieldData/FD/"+ss.str(),"maxSAR_idx", idx)==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
	}

	out_file.WriteAtrribute("/FieldData/FD","valid_cubes",m_Used);
	out_file.WriteAtrribute("/FieldData/FD","used_cubes",m_Valid);
	out_file.WriteAtrribute("/FieldData/FD","unused_cubes",m_Unused);
	out_file.WriteAtrribute("/FieldData/FD","air_cubes",m_AirVoxel);

	out_file.WriteAtrribute("/","proc_time", m_duration);
	return true;
}

bool SAR_Calculation::CalcFromHDF5(std::string h5_fn, std::string out_name, bool legacyHDF5)
{
	if (!ReadFromHDF5(h5_fn))
		return false;
	if (!CalcSAR())
		return false;
	return WriteToHDF5(out_name, legacyHDF5);
}


void SAR_Calculation::InitSAR()
{
	// delete any already existing results
	for (size_t i=0;i<m_SAR.size();++i)
		delete m_SAR.at(i);
	m_SAR.clear();

	DoAutoRange();

	m_posInOuput.Init("posInOutput", m_numLines);
	// this will fill m_posInOuput with all 0, 0 meaning not in the output data
	// a value > 0 gives the position+1 in the 3D output array
	unsigned int pos[3];
	unsigned int loc=0;
	for (unsigned int i=0; i<m_cellIndicies[0].size();++i)
	{
		pos[0] = m_cellIndicies[0].at(i);
		for (unsigned int j=0; j<m_cellIndicies[1].size();++j)
		{
			pos[1] = m_cellIndicies[1].at(j);
			for (unsigned int k=0; k<m_cellIndicies[2].size();++k)
			{
				pos[2] = m_cellIndicies[2].at(k);
				m_posInOuput({pos[0],pos[1],pos[2]}) = loc + 1;
				++loc;
			}
		}
	}

	unsigned int out_num_lines[3] = {(unsigned int)m_cellIndicies[0].size(),(unsigned int)m_cellIndicies[1].size(),(unsigned int)m_cellIndicies[2].size()};
	if ((m_record_cube_stats) && (m_freq.size()==1))
	{
		cout << "Enable Cube Statistics" << endl;
		m_cube_type.Init("CubeType", out_num_lines);
		m_cube_mass.Init("CubeMass", out_num_lines);
		m_cube_volume.Init("CubeVol", out_num_lines);
	}
	else if ((m_record_cube_stats) && (m_freq.size()>1))
		cerr << "Unable to store cube statistics for more than one frequency of interest!!" << endl;

	for (size_t n=0;n<m_freq.size();++n)
	{
		ArrayLib::ArrayIJK<float>* sar = new ArrayLib::ArrayIJK<float>("sar", out_num_lines);
		m_SAR.push_back(sar);
	}
	m_maxSAR.resize(m_freq.size(), 0);
	m_maxSAR_Idx.resize(m_freq.size());
}

void SAR_Calculation::DoAutoRange()
{
	if (m_autoRange<=0)
		return;
	if (m_DebugLevel>0)
		cout << "Calculate auto range with " << m_autoRange << "dB from loc max." << endl;

	std::array<unsigned int, 3> min_idx, max_idx;
	min_idx = {m_numLines[0],m_numLines[1],m_numLines[2]};
	max_idx = {0,0,0};
	unsigned int loc=0;
	float lim_fac = std::pow(10.0, m_autoRange / -10.0);
	std::vector<float> lim_pow;
	for (size_t n=0;n<m_freq.size();++n)
		lim_pow.push_back(lim_fac*m_local_cell_max_power_density.at(n));

	for (unsigned int i=0; i<m_numLines[0];++i)
		for (unsigned int j=0; j<m_numLines[1];++j)
			for (unsigned int k=0; k<m_numLines[2];++k)
			{
				for (size_t n=0;n<m_freq.size();++n)
				{
					if (m_local_cell_power_density.at(n)->data(loc) > lim_pow.at(n))
					{
						min_idx[0] = min(min_idx[0], i);
						min_idx[1] = min(min_idx[1], j);
						min_idx[2] = min(min_idx[2], k);
						max_idx[0] = max(max_idx[0], i+1);
						max_idx[1] = max(max_idx[1], j+1);
						max_idx[2] = max(max_idx[2], k+1);
					}
				}
				++loc;
			}
	for (int n=0;n<3;++n)
		SetSubRange(n, min_idx[n], max_idx[n]);
	// cerr << "DoAutoRange m_numLines" << m_numLines[0] << ", " << m_numLines[1] << ", " << m_numLines[2] << endl;
	// cerr << "DoAutoRange min" << min_idx[0] << ", " << min_idx[1] << ", " << min_idx[2] << endl;
	// cerr << "DoAutoRange max" << max_idx[0] << ", " << max_idx[1] << ", " << max_idx[2] << endl;
	}

bool SAR_Calculation::CalcLocalSAR()
{
	m_Valid=0;
	m_Used=0;
	m_Unused=0;
	m_AirVoxel=0;

	unsigned int N = m_cellIndicies[0].size() * m_cellIndicies[1].size() * m_cellIndicies[2].size();
	unsigned int loc;
	float dens;
	float lSAR;
	unsigned int pos[3];
	for (unsigned int i=0; i<m_cellIndicies[0].size();++i)
	{
		pos[0] = m_cellIndicies[0].at(i);
		for (unsigned int j=0; j<m_cellIndicies[1].size();++j)
		{
			pos[1] = m_cellIndicies[1].at(j);
			for (unsigned int k=0; k<m_cellIndicies[2].size();++k)
			{
				pos[2] = m_cellIndicies[2].at(k);
				dens = (*m_cell_density)(pos[0], pos[1], pos[2]);

				if (dens>0)
				{
					++m_Valid;
					for (size_t n=0;n<m_freq.size();++n)
					{
						lSAR = (*m_local_cell_power_density.at(n))(pos[0], pos[1], pos[2])/dens;
						(*m_SAR.at(n))(i,j,k) = lSAR;
						if (lSAR>m_maxSAR.at(n))
						{
							m_maxSAR.at(n) = lSAR;
							m_maxSAR_Idx.at(n) = {i,j,k};
						}
					}
				}
				else
				{
					++m_AirVoxel;
					// SAR(pos[0], pos[1], pos[2]) = 0;
				}
			}
		}
	}

	return true;
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
	double tol_mass = m_massTolerance*m_avg_mass;
	float min_size = 0;
	float max_size = 1e99;

	//iterate over cubical sizes to find fitting volume to mass
	while (mass_iterations<m_maxMassIterations)
	{
		// calculate cubical mass
		face_valid = GetCubicalMass(pos, box_size/2, start, stop, partial_start, partial_stop, mass, volume, bg_ratio, disabledFace);

		// check if found mass is valid
		mass_valid = abs(mass-m_avg_mass)<=tol_mass;
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
		if ((mass_iterations>0) && (face_valid==false) && (mass == old_mass) && (box_size>old_box_size) && (bg_ratio>=m_maxBGRatio))
		{
			// this is an invalid cube with a mass that cannot increase anymore
			// if (ignoreFaceValid)
			// 	return 0;
			return 2;
		}

		// record the lowest and largest box size as limits
		if (mass<m_avg_mass)
			min_size = max(min_size, box_size);
		else if (mass>m_avg_mass)
			max_size = min(max_size, box_size);


		// if no valid or finally invalid cube is found, calculate an alternative cube size
		if (mass_iterations==0)
		{
			// on first iteration, try a relative resize
			old_box_size=box_size;
			box_size*=pow(m_avg_mass/mass,1.0/3.0);
		}
		else if (mass == old_mass)
		{
			// if the mass happens to not have changed, try a value in the middle of known limits
			old_box_size=box_size;
			if (mass > m_avg_mass)
				box_size = 0.5*(min_size + box_size);
			else
				box_size = min(0.5*(max_size + box_size), 2.0*box_size); // plus make sure max_size is not still inf
		}
		else
		{
			// on later iterations, try a newton approach
			float new_box_size = box_size - (mass-m_avg_mass)/(mass-old_mass)*(box_size-old_box_size);
			if (new_box_size<min_size) // new box size below lower limit
				new_box_size = 0.5*(min_size + box_size);
			else if (new_box_size>max_size) // new box size above upper limit
				new_box_size = 0.5*(max_size + box_size);
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
		cerr << "SAR_Calculation::GetCubicalMass: critical error: invalid averaging box size: " << box_size << " !! EXIT" << endl;
		exit(-1);
	}
	bool face_valid=true;
	double n_box_size;
	for (int n=0;n<3;++n)
	{
		n_box_size = box_size;
		// check start position
		start[n]=pos[n];
		partial_start[n]=1;
		float dist=m_cellWidth[n][pos[n]]/2;
		if ((disabledFace==2*n) || (disabledFace==2*n+1))
			n_box_size = 2*box_size - m_cellWidth[n][pos[n]]/2;

		if (disabledFace!=2*n)
		{
			while (dist<n_box_size)
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
				partial_start[n]=1-(dist-n_box_size)/m_cellWidth[n][start[n]];
			}
			if ((partial_start[n]!=-1) && (pos[n]==start[n]))
				partial_start[n]=2*n_box_size/m_cellWidth[n][start[n]];
		}

		// check stop position
		stop[n]=pos[n];
		partial_stop[n]=1;
		dist=m_cellWidth[n][pos[n]]/2;
		if (disabledFace!=2*n+1)
		{
			while (dist<n_box_size)
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
				partial_stop[n]=1-(dist-n_box_size)/m_cellWidth[n][stop[n]];
			}
			if ((partial_stop[n]!=-1) && (pos[n]==stop[n]))
				partial_stop[n]=2*n_box_size/m_cellWidth[n][stop[n]];
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
	double cell_vol=0;
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

				cell_vol = (*m_cell_volume)(f_pos[0], f_pos[1], f_pos[2])*weight[0]*weight[1]*weight[2];

				mass += (*m_cell_density)(f_pos[0], f_pos[1], f_pos[2])*cell_vol;
				volume += cell_vol;

				if ((*m_cell_density)(f_pos[0], f_pos[1], f_pos[2])==0)
					bg_volume += cell_vol;
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
		face_valid = face_valid && face_intersect[n];

	bg_ratio = bg_volume/volume;

	return face_valid;
}

void SAR_Calculation::CalcCubicalSAR(double* vx_sar, unsigned int pos[3], unsigned int start[3], unsigned int stop[3],
									  float partial_start[3], float partial_stop[3],
									  ArrayLib::ArrayIJK<bool> &Vx_Used, ArrayLib::ArrayIJK<bool> &Vx_Valid,
									  bool assignUsed)
{
	double power_mass=0;
	double mass=0;
	double weight[3];
	unsigned int f_pos[3];
	unsigned int loc;
	size_t N=m_SAR.size();
	for (size_t n=0;n<N;++n)
		vx_sar[n] = 0;
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

				loc = m_cell_density->linearIndex({f_pos[0], f_pos[1], f_pos[2]});
				if (m_cell_density->data(loc)>=0)
				{
					mass += m_cell_density->data(loc)*m_cell_volume->data(loc)*weight[0]*weight[1]*weight[2];
					for (size_t n=0;n<N;++n)
						vx_sar[n] += m_local_cell_power_density.at(n)->data(loc)*m_cell_volume->data(loc)*weight[0]*weight[1]*weight[2];
				}
			}
		}
	}
	for (size_t n=0;n<N;++n)
		vx_sar[n] /= mass;
	// float vx_SAR = power_mass/mass;
	if (assignUsed==false)
		return;


	unsigned int out_loc = m_posInOuput(pos[0], pos[1], pos[2]);
	if (out_loc==0)
		// this should not happen!
		throw;

	for (size_t n=0;n<N;++n)
		m_SAR.at(n)->data(out_loc-1)=vx_sar[n];

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
					loc = m_cell_density->linearIndex({f_pos[0], f_pos[1], f_pos[2]});
					out_loc = m_posInOuput.data(loc);
					if (out_loc==0)
						continue;
					--out_loc;
					if (!Vx_Valid.data(out_loc) && (m_cell_density->data(loc)>0))
					{
						Vx_Used.data(out_loc)=true;
						for (size_t n=0;n<N;++n)
							m_SAR.at(n)->data(out_loc)=max(m_SAR.at(n)->data(out_loc), (float)vx_sar[n]);
					}
				}
			}
		}
	}
	return;
}


bool SAR_Calculation::CalcAveragedSAR()
{
	unsigned int pos[3];

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

	unsigned int out_num_lines[3] = {(unsigned int)m_cellIndicies[0].size(),(unsigned int)m_cellIndicies[1].size(),(unsigned int)m_cellIndicies[2].size()};
	ArrayLib::ArrayIJK<bool> Vx_Used("vx_used", out_num_lines);
	bool* vx_used_data = Vx_Used.data();
	ArrayLib::ArrayIJK<bool> Vx_Valid("vx_used", out_num_lines);
	bool* vx_valid_data = Vx_Valid.data();
	unsigned char* vx_cube_type = m_cube_type.data();
	float* cube_mass =  m_cube_mass.data();
	float* cube_vol = m_cube_volume.data();
	float* cell_density = m_cell_density->data();
	unsigned int loc = 0;
	unsigned int out_loc = 0;
	unsigned int prog_N = max(1,(int)(Vx_Valid.size()/100));
	double* vx_sar = new double[m_SAR.size()];
	if (m_progress)
		std::cout << "Step-1: ";

	for (unsigned int i=0; i<m_cellIndicies[0].size();++i)
	{
		pos[0] = m_cellIndicies[0].at(i);
		for (unsigned int j=0; j<m_cellIndicies[1].size();++j)
		{
			pos[1] = m_cellIndicies[1].at(j);
			for (unsigned int k=0; k<m_cellIndicies[2].size();++k)
			{
				pos[2] = m_cellIndicies[2].at(k);
				loc = m_cell_density->linearIndex({pos[0],pos[1],pos[2]});
				out_loc = Vx_Valid.linearIndex({i,j,k});
				if (cell_density[loc]==0)
				{
					// SAR(i,j,k) = 0;
					++m_AirVoxel;
					// ++loc;
					if (m_progress && (loc%prog_N == 0))
					{
						std::cout << ".";
						std::cout.flush();
					}
					continue;
				}

				// guess an initial box size and find a fitting cube
				EC = FindFittingCubicalMass(pos, pow(m_avg_mass/cell_density[loc],1.0/3.0)/2, start, stop,
											partial_start, partial_stop, total_mass, voxel_volume, bg_ratio, -1, m_IgnoreFaceValid);

				if (EC==0)
				{
					vx_valid_data[out_loc]=true;
					vx_used_data[out_loc]=true;
					if (vx_cube_type!=NULL)
						vx_cube_type[out_loc]=7;
					if (cube_mass!=NULL)
						cube_mass[out_loc]=total_mass;
					if (cube_vol!=NULL)
						cube_vol[out_loc]=voxel_volume;
					++m_Valid;
					CalcCubicalSAR(vx_sar, pos, start, stop, partial_start, partial_stop, Vx_Used, Vx_Valid, true);
					for (int n=0;n<m_maxSAR.size();++n)
						if (vx_sar[n]>m_maxSAR.at(n))
						{
							m_maxSAR.at(n) = vx_sar[n];
							m_maxSAR_Idx.at(n) = {i,j,k};
						}
				}
				else if (EC==1)
					++cnt_case1;
				else if (EC==2)
					++cnt_case2;
				else if (EC==-1)
					++cnt_NoConvergence;
				else
					cerr << "other EC" << EC << endl;
				// ++loc;
				if (m_progress && (out_loc%prog_N == 0))
				{
					std::cout << ".";
					std::cout.flush();
				}
			}
		}
	}
	delete[] vx_sar;
	if (m_progress)
		std::cout << " Done" << std::endl;

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

	ArrayLib::ArrayIJK<float> dummy_SAR;
	// count all used and unused etc. + special handling of unused voxels!!
	m_Used=0;
	m_Unused=0;
	loc = 0;
	if (m_progress)
		std::cout << "Step-2: ";

	unsigned int N = m_SAR.size();
	double* vx_sar_unused[6];
	for (int n=0;n<6;++n)
		vx_sar_unused[n] = new double[N];

	for (unsigned int i=0; i<m_cellIndicies[0].size();++i)
	{
		pos[0] = m_cellIndicies[0].at(i);
		for (unsigned int j=0; j<m_cellIndicies[1].size();++j)
		{
			pos[1] = m_cellIndicies[1].at(j);
			for (unsigned int k=0; k<m_cellIndicies[2].size();++k)
			{
				pos[2] = m_cellIndicies[2].at(k);
				loc = m_cell_density->linearIndex({pos[0],pos[1],pos[2]});
				out_loc = Vx_Valid.linearIndex({i,j,k});
				if (!vx_valid_data[out_loc] && vx_used_data[out_loc])
				{
					++m_Used;
					if (vx_cube_type!=NULL)
						vx_cube_type[out_loc] = 8;
				}
				if ((cell_density[loc]>0) && !vx_valid_data[out_loc] && !vx_used_data[out_loc])
				{
					++m_Unused;

					for (int n=0;n<N;++n)
						m_SAR.at(n)->data(out_loc) = 0;
					double unused_volumes[6];
					double unused_mass[6];

					double min_Vol=FLT_MAX;

					// special handling of unused voxels:
					for (int n=0;n<6;++n)
					{
						EC = FindFittingCubicalMass(pos, pow(m_avg_mass/cell_density[loc],1.0/3.0)/2, start, stop,
													partial_start, partial_stop, unused_mass[n], unused_volumes[n], bg_ratio, n, true);
						if (EC==0)
						{
							// cerr << "calc unused SAR " << n << ": " << pos[n/2] << " " << start[n/2] << "; " << partial_start[n/2] << " " << stop[n/2] << " " << partial_stop[n/2] << endl;
							CalcCubicalSAR(vx_sar_unused[n], pos, start, stop, partial_start, partial_stop, Vx_Used, Vx_Valid, false);
							min_Vol = min(min_Vol,unused_volumes[n]);
						}
						else if (EC==2)
						{
							// target mass was not reached!
							for (int m=0;m<N;++m)
								vx_sar_unused[n][m] = 0;
						}
						else
						{
							// this should not happen
							cerr << EC << " SAR_Calculation::CalcAveragedSAR: Error handling unused voxels, can't find fitting cubical averaging volume' " << endl;
							for (int m=0;m<N;++m)
								vx_sar_unused[n][m] = 0;
						}
					}
					for (int n=0;n<6;++n)
					{
						if (unused_volumes[n]<=m_UnusedRelativeVolLimit*min_Vol)
							for (int m=0;m<N;++m)
							{
								if (vx_sar_unused[n][m]>m_SAR.at(m)->data(out_loc))
								{
									if (vx_sar_unused[n][m]>m_maxSAR.at(m))
									{
										m_maxSAR.at(m) = vx_sar_unused[n][m];
										m_maxSAR_Idx.at(m) = {i,j,k};
									}
									m_SAR.at(m)->data(out_loc) = vx_sar_unused[n][m];
									if (vx_cube_type!=NULL)
										vx_cube_type[out_loc] = n+1;
									if (cube_mass!=NULL)
										cube_mass[out_loc]=unused_mass[n];
									if (cube_vol!=NULL)
										cube_vol[out_loc]=unused_volumes[n];
								}
							}
					}
				}
				if (m_progress && (out_loc%prog_N == 0))
				{
					std::cout << ".";
					std::cout.flush();
				}
			}
		}
	}
	for (int n=0;n<6;++n)
		delete[] vx_sar_unused[n];

	if (m_progress)
		std::cout << " Done" << std::endl;

	if (m_Valid+m_Used+m_Unused+m_AirVoxel!=Vx_Valid.size())
	{
		cerr << "SAR_Calculation::CalcAveragedSAR: critical error, mismatch in voxel status count... EXIT" << endl;
		exit(1);
	}

	if (m_DebugLevel>0)
		cerr << "SAR_Calculation::CalcAveragedSAR: Stats: Valid=" << m_Valid << " Used=" << m_Used << " Unused=" << m_Unused << " Air-Voxel=" << m_AirVoxel << endl;

	return true;
}
