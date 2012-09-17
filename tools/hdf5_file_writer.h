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

class HDF5_File_Writer
{
public:
	HDF5_File_Writer(std::string filename);
	~HDF5_File_Writer();

	bool WriteRectMesh(unsigned int const* numLines, double const* const* discLines, int MeshType=0, double scaling=1);
	bool WriteRectMesh(unsigned int const* numLines, float const* const* discLines, int MeshType=0, float scaling=1);

	bool WriteScalarField(std::string dataSetName, float const* const* const* field, size_t datasize[3]);
	bool WriteScalarField(std::string dataSetName, double const* const* const* field, size_t datasize[3]);

	bool WriteScalarField(std::string dataSetName, std::complex<float> const* const* const* field, size_t datasize[3]);
	bool WriteScalarField(std::string dataSetName, std::complex<double> const* const* const* field, size_t datasize[3]);

	bool WriteVectorField(std::string dataSetName, float const* const* const* const* field, size_t datasize[3]);
	bool WriteVectorField(std::string dataSetName, double const* const* const* const* field, size_t datasize[3]);

	bool WriteVectorField(std::string dataSetName, std::complex<float> const* const* const* const* field, size_t datasize[3]);
	bool WriteVectorField(std::string dataSetName, std::complex<double> const* const* const* const* field, size_t datasize[3]);

	bool WriteData(std::string dataSetName, float const* field_buf, size_t dim, size_t* datasize);
	bool WriteData(std::string dataSetName, double const* field_buf, size_t dim, size_t* datasize);

	bool WriteAtrribute(std::string locName, std::string attr_name, void const* value, hsize_t size, hid_t mem_type);
	bool WriteAtrribute(std::string locName, std::string attr_name, float const* value, hsize_t size);
	bool WriteAtrribute(std::string locName, std::string attr_name, double const* value, hsize_t size);
	bool WriteAtrribute(std::string locName, std::string attr_name, std::vector<float> values);
	bool WriteAtrribute(std::string locName, std::string attr_name, std::vector<double> values);
	bool WriteAtrribute(std::string locName, std::string attr_name, float value);
	bool WriteAtrribute(std::string locName, std::string attr_name, double value);

	void SetCurrentGroup(std::string group, bool createGrp=true);

protected:
	std::string m_filename;
	std::string m_Group;

	hid_t OpenGroup(hid_t hdf5_file, std::string group);
	bool WriteData(std::string dataSetName, hid_t mem_type, void const* field_buf, size_t dim, size_t* datasize);
};

#endif // HDF5_FILE_WRITER_H
