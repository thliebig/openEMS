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

#ifndef HDF5_FILE_WRITER_H
#define HDF5_FILE_WRITER_H

#include <string>
#include <vector>
#include <complex>
#include <hdf5.h>

#include "arraylib/array_nijk.h"
#include "arraylib/array_ijk.h"

#define OPENEMS_HDF5_VERSION 0.3

class HDF5_File_Writer
{
public:
	HDF5_File_Writer(std::string filename);
	~HDF5_File_Writer();

	bool WriteRectMesh(unsigned int const* numLines, double const* const* discLines, int MeshType=0, double scaling=1, std::string s_mesh_grp="/Mesh");
	bool WriteRectMesh(unsigned int const* numLines, float const* const* discLines, int MeshType=0, double scaling=1, std::string s_mesh_grp="/Mesh");

	template <typename T>
	bool WriteScalarField(std::string dataSetName, ArrayLib::ArrayIJK<T> &data, bool legacy_fmt=false);

	template <typename T>
	bool WriteVectorField(std::string dataSetName, ArrayLib::ArrayNIJK<T> &data, bool legacy_fmt=false);


	bool WriteData(std::string dataSetName, float const* field_buf, size_t dim, size_t* datasize, std::string d_order="");
	bool WriteData(std::string dataSetName, double const* field_buf, size_t dim, size_t* datasize, std::string d_order="");
	bool WriteData(std::string dataSetName, std::complex<float> const* field_buf, size_t dim, size_t* datasize, std::string d_order="");
	bool WriteData(std::string dataSetName, std::complex<double> const* field_buf, size_t dim, size_t* datasize, std::string d_order="");

	bool WriteAtrribute(std::string locName, std::string attr_name, float const* value, hsize_t size);
	bool WriteAtrribute(std::string locName, std::string attr_name, double const* value, hsize_t size);
	bool WriteAtrribute(std::string locName, std::string attr_name, std::vector<float> values);
	bool WriteAtrribute(std::string locName, std::string attr_name, std::vector<double> values);
	bool WriteAtrribute(std::string locName, std::string attr_name, float value);
	bool WriteAtrribute(std::string locName, std::string attr_name, double value);
	bool WriteAtrribute(std::string locName, std::string attr_name, int value);
	bool WriteAtrribute(std::string locName, std::string attr_name, unsigned int value);
	bool WriteAtrribute(std::string locName, std::string attr_name, std::string value);
	bool WriteAtrribute(std::string locName, std::string attr_name, bool value);

	void SetCurrentGroup(std::string group, bool createGrp=true);

	bool WriteData(std::string dataSetName, hid_t mem_type, void const* field_buf, size_t dim, size_t* datasize, std::string d_order="");
	bool WriteAtrribute(std::string locName, std::string attr_name, void const* value, hsize_t size, hid_t mem_type);

protected:
	std::string m_filename;
	std::string m_Group;

	hid_t OpenGroup(hid_t hdf5_file, std::string group);
	bool WriteData(hid_t group, std::string dataSetName, hid_t mem_type, void const* field_buf, size_t dim, size_t* datasize, std::string d_order="");
	bool WriteAtrribute(hid_t loc, std::string attr_name, void const* value, hsize_t size, hid_t mem_type);
};

#endif // HDF5_FILE_WRITER_H
