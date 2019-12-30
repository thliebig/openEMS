/*
*	Copyright (C) 2011 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef HDF5_FILE_READER_H
#define HDF5_FILE_READER_H

#include <string>
#include <vector>
#include <complex>
#include <hdf5.h>
#define _USE_MATH_DEFINES

class HDF5_File_Reader
{
public:
	HDF5_File_Reader(std::string filename);
	virtual ~HDF5_File_Reader();

	bool ReadMesh(float** lines, unsigned int* numLines, int &meshType);

	//! Get the number of timesteps stored at /FieldData/TD/<NUMBER_OF_TS>
	unsigned int GetNumTimeSteps();
	bool ReadTimeSteps(std::vector<unsigned int> &timestep, std::vector<std::string> &names);

	/*!
	  Get time-domain data stored at /FieldData/TD/<NUMBER_OF_TS>
	  \param[in]  idx	time step index to extract
	  \param[out] time	time attribute for the given timestep
	  \param[out] data_size data size found
	  \return field data found in given timestep, caller must delete array, returns NULL if timestep was not found
	  */
	float**** GetTDVectorData(size_t idx, float &time, unsigned int data_size[4]);

	unsigned int GetNumFrequencies();
	bool ReadFrequencies(std::vector<float> &frequencies);
	bool ReadFrequencies(std::vector<double> &frequencies);

	/*!
	  Get frequency-domain data stored at "/FieldData/FD/f<idx>_real" and "/FieldData/FD/f<idx>_imag"
	  \param[in]  idx	frequency index to extract
	  \param[out] data_size data size found
	  \return complex field data found for the given frequency, caller must delete array, returns NULL if frequency was not found
	  */
	std::complex<float>**** GetFDVectorData(size_t idx, unsigned int data_size[4]);

	/*!
	  Calculate
	  */
	bool CalcFDVectorData(std::vector<float> &frequencies, std::vector<std::complex<float>****> &FD_data, unsigned int data_size[4]);

	bool ReadAttribute(std::string grp_name, std::string attr_name, std::vector<double> &attr_values);
	bool ReadAttribute(std::string grp_name, std::string attr_name, std::vector<float> &attr_values);

	bool IsValid();

protected:
	std::string m_filename;

	bool ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, double* &data);
	bool ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, float* &data);

	bool OpenGroup(hid_t &file, hid_t &group, std::string groupName);
};

#endif // HDF5_FILE_READER_H
