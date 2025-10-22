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
#include "arraylib/array_nijk.h"
#include "arraylib/array_ijk.h"

#define _USE_MATH_DEFINES

class HDF5_File_Reader
{
public:
	HDF5_File_Reader(std::string filename);
	virtual ~HDF5_File_Reader();

	//! Check if a legacy openEMS file format was detected
	bool GetLegacyFormat() const {return m_legacy_fmt;}

	//! Get the openEMS file version as defined in the hdf5 file, 0 if not found
	double GetFileVersion() const {return m_file_ver;}

	bool ReadMesh(double** lines, unsigned int* numLines, int &meshType, std::string s_mesh_grp="/Mesh");
	bool ReadMesh(float** lines, unsigned int* numLines, int &meshType, std::string s_mesh_grp="/Mesh");

	//! Get the number of timesteps stored at /FieldData/TD/<NUMBER_OF_TS>
	bool HasTDData() {return CheckDataSetExists("/FieldData/TD");}
	unsigned int GetNumTimeSteps();
	bool ReadTimeSteps(std::vector<unsigned int> &timestep, std::vector<std::string> &names);

	/*!
	  Get time-domain data stored at /FieldData/TD/<NUMBER_OF_TS>
	  \param[in]  idx	time step index to extract
	  \param[out] time	time attribute for the given timestep
	  \param[out] data  reference to 4d array for data storage
	  \return true if successful
	  */
	bool GetTDVectorData(size_t idx, float &time, ArrayLib::ArrayNIJK<float> &data);

	bool HasFDData() {return CheckDataSetExists("/FieldData/FD");}
	unsigned int GetNumFrequencies();
	bool ReadFrequencies(std::vector<float> &frequencies);
	bool ReadFrequencies(std::vector<double> &frequencies);

	/*!
	  Get frequency-domain data stored at "/FieldData/FD/f<idx>_real" and "/FieldData/FD/f<idx>_imag"
	  \param[in]  idx  frequency index to extract
	  \param[out] data reference to a complex 4D arrays
	  \return true if successful
	  */
	bool GetFDVectorData(size_t idx, ArrayLib::ArrayNIJK<std::complex<float>> &data);

	template <typename T>
	bool ReadVectorDataSet(std::string ds_name, ArrayLib::ArrayNIJK<T> &data, bool type_check=true);

	/*!
	  Calculate frequency-domain vector data based on an array of time domain data using a DFT.
	  \param[in]  frequencies  list of frequencies to calculate
	  \param[out] FD_data      reference to a list of 4D arrays. Caller must cleanup the array data
	  \return true if successful
	  */
	bool CalcFDVectorData(std::vector<float> &frequencies, std::vector<ArrayLib::ArrayNIJK<std::complex<float>>*> &FD_data);

	bool ReadAttribute(std::string grp_name, std::string attr_name, std::vector<double> &attr_values);
	bool ReadAttribute(std::string grp_name, std::string attr_name, std::vector<float> &attr_values);

	template <typename T>
	bool ReadAttribute(std::string grp_name, std::string attr_name, T &value, bool silent=false);

	template <typename T>
	bool ReadScalarDataSet(std::string ds_name, ArrayLib::ArrayIJK<T> &data, bool type_check=true);

	bool GetDataSetNameByIndex(std::string grp_name, unsigned int idx, std::string &name);

	bool CheckDataSetExists(std::string ds_name);
	bool IsValid();

	template <typename T>
	static hid_t GetH5Type();
	template <typename T>
	static bool CheckH5Type(hid_t type);


protected:
	std::string m_filename;
	hid_t m_hdf5_file;
	hid_t getFile();
	bool open();
	void close();
	double m_file_ver = 0;
	bool m_legacy_fmt = false;

	bool ReadDataSet(std::string ds_name, hid_t type, void* buffer);

	bool GetDataSetInfo(std::string ds_name, hid_t &type, hsize_t &nDim, hsize_t* &dims);

	template <typename T>
	bool ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, T* &data);

	bool ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, float* &f_data, double* &d_data);
	bool ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, double* &data);
	bool ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, float* &data);

	bool OpenGroup(hid_t &file, hid_t &group, std::string groupName);

	bool GetDataSetNameByIndex(hid_t &group, unsigned int idx, std::string &name);
};

#endif // HDF5_FILE_READER_H
