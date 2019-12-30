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

using namespace std;

#include "hdf5_file_writer.h"
#include <boost/algorithm/string.hpp>
#include <hdf5.h>

#include <sstream>
#include <iostream>
#include <iomanip>

HDF5_File_Writer::HDF5_File_Writer(string filename)
{
	m_filename = filename;
	m_Group = "/";
	hid_t hdf5_file = H5Fcreate(m_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (hdf5_file<0)
	{
		cerr << "HDF5_File_Writer::HDF5_File_Writer: Error, creating the given file """ << m_filename << """ failed" << endl;
	}
	H5Fclose(hdf5_file);
}

HDF5_File_Writer::~HDF5_File_Writer()
{
}

hid_t HDF5_File_Writer::OpenGroup(hid_t hdf5_file, string group)
{
	if (hdf5_file<0)
	{
		cerr << "HDF5_File_Writer::CreateGroup: Error, invalid file id" << endl;
		return -1;
	}

	vector<string> results;
	boost::split(results, group, boost::is_any_of("/"));

	hid_t grp=H5Gopen(hdf5_file,"/", H5P_DEFAULT);
	if (grp<0)
	{
		cerr << "HDF5_File_Writer::OpenGroup: Error, opening root group " << endl;
		return -1;
	}

	for (size_t n=0;n<results.size();++n)
	{
		hid_t old_grp = grp;
		if (!results.at(n).empty())
		{
			if (H5Lexists(grp, results.at(n).c_str(), H5P_DEFAULT))
			{
			  grp = H5Gopen(grp, results.at(n).c_str(), H5P_DEFAULT);
				H5Gclose(old_grp);
				if (grp<0)
				{
					cerr << "HDF5_File_Writer::OpenGroup: Error, failed to open existing group" << endl;
					return -1;
				}
			}
			else
			{
			  grp = H5Gcreate(grp,results.at(n).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				H5Gclose(old_grp);
				if (grp<0)
				{
					cerr << "HDF5_File_Writer::OpenGroup: Error, creating group """ << group << """ failed" << endl;
					return -1;
				}
			}
		}
	}
	return grp;
}

void HDF5_File_Writer::SetCurrentGroup(std::string group, bool createGrp)
{
	m_Group = group;
	if (createGrp==false)
		return;

	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	if (hdf5_file<0)
	{
		cerr << "HDF5_File_Writer::SetCurrentGroup: Error, opening the given file """ << m_filename << """ failed" << endl;
		return;
	}
	hid_t hdf5_group = OpenGroup(hdf5_file, m_Group);
	if (hdf5_group<0)
		cerr << "HDF5_File_Writer::WriteData: Error opening group" << endl;
	H5Gclose(hdf5_group);
	H5Fclose(hdf5_file);
}

bool HDF5_File_Writer::WriteRectMesh(unsigned int const* numLines, double const* const* discLines, int MeshType, double scaling)
{
	float* array[3];
	for (int n=0; n<3; ++n)
	{
		array[n] = new float[numLines[n]];
		for (unsigned int i=0; i<numLines[n]; ++i)
			array[n][i]=discLines[n][i];

	}
	bool success = WriteRectMesh(numLines,array,MeshType,scaling);
	for (int n=0; n<3; ++n)
		delete[] array[n];
	return success;
}

bool HDF5_File_Writer::WriteRectMesh(unsigned int const* numLines, float const* const* discLines, int MeshType, float scaling)
{
	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	if (hdf5_file<0)
	{
		cerr << "HDF5_File_Writer::WriteRectMesh: Error, opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	if (H5Lexists(hdf5_file, "/Mesh", H5P_DEFAULT))
	{
		cerr << "HDF5_File_Writer::WriteRectMesh: Error, group ""/Mesh"" already exists" << endl;
		H5Fclose(hdf5_file);
		return false;
	}

	hid_t mesh_grp = H5Gcreate(hdf5_file,"/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (mesh_grp<0)
	{
		cerr << "HDF5_File_Writer::WriteRectMesh: Error, creating group ""/Mesh"" failed" << endl;
		H5Fclose(hdf5_file);
		return false;
	}

	string names[] = {"x","y","z"};
	if (MeshType==1)
	{
		names[0]="rho";
		names[1]="alpha";
	}
	if (MeshType==2)
	{
		names[0]="r";
		names[1]="theta";
		names[2]="phi";
	}

	for (int n=0; n<3; ++n)
	{
		hsize_t dims[1]={numLines[n]};
		hid_t space = H5Screate_simple(1, dims, NULL);
		hid_t dataset = H5Dcreate(mesh_grp, names[n].c_str(), H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		float* array = new float[numLines[n]];
		for (unsigned int i=0; i<numLines[n]; ++i)
		{
			if ((MeshType==1) && (n==1)) //check for alpha-direction
				array[i] = discLines[n][i];
			else if ((MeshType==2) && (n>0)) //check for theta/phi-direction
				array[i] = discLines[n][i];
			else
				array[i] = discLines[n][i] * scaling;
		}
		if (H5Dwrite(dataset, H5T_NATIVE_FLOAT, space, H5P_DEFAULT, H5P_DEFAULT, array))
		{
			cerr << "HDF5_File_Writer::WriteRectMesh: Error, writing to dataset failed" << endl;
			delete[] array;
			H5Dclose(dataset);
			H5Sclose(space);
			H5Gclose(mesh_grp);
			H5Fclose(hdf5_file);
			return false;
		}
		delete[] array;
		H5Dclose(dataset);
		H5Sclose(space);
	}
	H5Gclose(mesh_grp);
	H5Fclose(hdf5_file);
	return true;
}

bool HDF5_File_Writer::WriteScalarField(std::string dataSetName, float const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2];
	size_t n_size[3]={datasize[2],datasize[1],datasize[0]};
	float* buffer = new float[size];
	for (size_t k=0;k<datasize[2];++k)
		for (size_t j=0;j<datasize[1];++j)
			for (size_t i=0;i<datasize[0];++i)
			{
				buffer[pos++]=field[i][j][k];
			}
	bool success = WriteData(dataSetName,buffer,3,n_size);
	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteScalarField(std::string dataSetName, double const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2];
	size_t n_size[3]={datasize[2],datasize[1],datasize[0]};
	double* buffer = new double[size];
	for (size_t k=0;k<datasize[2];++k)
		for (size_t j=0;j<datasize[1];++j)
			for (size_t i=0;i<datasize[0];++i)
			{
				buffer[pos++]=field[i][j][k];
			}
	bool success = WriteData(dataSetName,buffer,3,n_size);
	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteScalarField(std::string dataSetName, std::complex<float> const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2];
	size_t n_size[3]={datasize[2],datasize[1],datasize[0]};
	float* buffer = new float[size];
	for (size_t k=0;k<datasize[2];++k)
		for (size_t j=0;j<datasize[1];++j)
			for (size_t i=0;i<datasize[0];++i)
			{
				buffer[pos++]=real(field[i][j][k]);
			}
	bool success = WriteData(dataSetName + "_real",buffer,3,n_size);

	pos = 0;
	for (size_t k=0;k<datasize[2];++k)
		for (size_t j=0;j<datasize[1];++j)
			for (size_t i=0;i<datasize[0];++i)
			{
				buffer[pos++]=imag(field[i][j][k]);
			}
	success &= WriteData(dataSetName + "_imag",buffer,3,n_size);

	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteScalarField(std::string dataSetName, std::complex<double> const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2];
	size_t n_size[3]={datasize[2],datasize[1],datasize[0]};
	double* buffer = new double[size];
	for (size_t k=0;k<datasize[2];++k)
		for (size_t j=0;j<datasize[1];++j)
			for (size_t i=0;i<datasize[0];++i)
			{
				buffer[pos++]=real(field[i][j][k]);
			}
	bool success = WriteData(dataSetName + "_real",buffer,3,n_size);

	pos = 0;
	for (size_t k=0;k<datasize[2];++k)
		for (size_t j=0;j<datasize[1];++j)
			for (size_t i=0;i<datasize[0];++i)
			{
				buffer[pos++]=imag(field[i][j][k]);
			}
	success &= WriteData(dataSetName + "_imag",buffer,3,n_size);

	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteVectorField(std::string dataSetName, float const* const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2]*3;
	float* buffer = new float[size];
	for (int n=0;n<3;++n)
		for (size_t k=0;k<datasize[2];++k)
			for (size_t j=0;j<datasize[1];++j)
				for (size_t i=0;i<datasize[0];++i)
				{
					buffer[pos++]=field[n][i][j][k];
				}
	size_t n_size[4]={3,datasize[2],datasize[1],datasize[0]};
	bool success = WriteData(dataSetName,buffer,4,n_size);
	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteVectorField(std::string dataSetName, double const* const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2]*3;
	double* buffer = new double[size];
	for (int n=0;n<3;++n)
		for (size_t k=0;k<datasize[2];++k)
			for (size_t j=0;j<datasize[1];++j)
				for (size_t i=0;i<datasize[0];++i)
				{
					buffer[pos++]=field[n][i][j][k];
				}
	size_t n_size[4]={3,datasize[2],datasize[1],datasize[0]};
	bool success = WriteData(dataSetName,buffer,4,n_size);
	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteVectorField(std::string dataSetName, std::complex<float> const* const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2]*3;
	size_t n_size[4]={3,datasize[2],datasize[1],datasize[0]};
	float* buffer = new float[size];
	for (int n=0;n<3;++n)
		for (size_t k=0;k<datasize[2];++k)
			for (size_t j=0;j<datasize[1];++j)
				for (size_t i=0;i<datasize[0];++i)
				{
					buffer[pos++]=real(field[n][i][j][k]);
				}
	bool success = WriteData(dataSetName + "_real",buffer,4,n_size);

	pos = 0;
	for (int n=0;n<3;++n)
		for (size_t k=0;k<datasize[2];++k)
			for (size_t j=0;j<datasize[1];++j)
				for (size_t i=0;i<datasize[0];++i)
				{
					buffer[pos++]=imag(field[n][i][j][k]);
				}
	success &= WriteData(dataSetName + "_imag",buffer,4,n_size);

	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteVectorField(std::string dataSetName, std::complex<double> const* const* const* const* field, size_t datasize[3])
{
	size_t pos = 0;
	size_t size = datasize[0]*datasize[1]*datasize[2]*3;
	size_t n_size[4]={3,datasize[2],datasize[1],datasize[0]};
	double* buffer = new double[size];
	for (int n=0;n<3;++n)
		for (size_t k=0;k<datasize[2];++k)
			for (size_t j=0;j<datasize[1];++j)
				for (size_t i=0;i<datasize[0];++i)
				{
					buffer[pos++]=real(field[n][i][j][k]);
				}
	bool success = WriteData(dataSetName + "_real",buffer,4,n_size);

	pos = 0;
	for (int n=0;n<3;++n)
		for (size_t k=0;k<datasize[2];++k)
			for (size_t j=0;j<datasize[1];++j)
				for (size_t i=0;i<datasize[0];++i)
				{
					buffer[pos++]=imag(field[n][i][j][k]);
				}
	success &= WriteData(dataSetName + "_imag",buffer,4,n_size);

	delete[] buffer;
	return success;
}

bool HDF5_File_Writer::WriteData(std::string dataSetName, float const* field_buf, size_t dim, size_t* datasize)
{
	return WriteData(dataSetName, H5T_NATIVE_FLOAT, field_buf,dim, datasize);
}

bool HDF5_File_Writer::WriteData(std::string dataSetName, double const* field_buf, size_t dim, size_t* datasize)
{
	return WriteData(dataSetName, H5T_NATIVE_DOUBLE, field_buf,dim, datasize);
}

bool HDF5_File_Writer::WriteData(std::string dataSetName,  hid_t mem_type, void const* field_buf, size_t dim, size_t* datasize)
{
	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	if (hdf5_file<0)
	{
		cerr << "HDF5_File_Writer::WriteData: Error, opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hid_t group = OpenGroup(hdf5_file,m_Group);
	if (group<0)
	{
		cerr << "HDF5_File_Writer::WriteData: Error opening group" << endl;
		H5Fclose(hdf5_file);
		return false;
	}

	hsize_t* dims = new hsize_t[dim];
	for (size_t n=0;n<dim;++n)
		dims[n]=datasize[n];
	hid_t space = H5Screate_simple(dim, dims, NULL);
	hid_t dataset = H5Dcreate(group, dataSetName.c_str(), mem_type, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (H5Dwrite(dataset, mem_type, space, H5P_DEFAULT, H5P_DEFAULT, field_buf))
	{
		cerr << "HDF5_File_Writer::WriteData: Error, writing to dataset failed" << endl;
		H5Dclose(dataset);
		H5Sclose(space);
		H5Gclose(group);
		H5Fclose(hdf5_file);
		return false;
	}
	H5Dclose(dataset);
	H5Sclose(space);
	H5Gclose(group);
	H5Fclose(hdf5_file);
	delete[] dims; dims=NULL;
	return true;
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, void const* value, hsize_t size, hid_t mem_type)
{
	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT );
	if (hdf5_file<0)
	{
		cerr << "HDF5_File_Writer::WriteAtrribute: Error, opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	if (H5Lexists(hdf5_file, locName.c_str(), H5P_DEFAULT)<0)
	{
		cerr << "HDF5_File_Writer::WriteAtrribute: Error, failed to find location: """ << locName << """" << endl;
		H5Fclose(hdf5_file);
		return false;
	}
	hid_t loc = H5Oopen(hdf5_file, locName.c_str(), H5P_DEFAULT);
	if (loc<0)
	{
		cerr << "HDF5_File_Writer::WriteAtrribute: Error, failed to open location: """ << locName << """" << endl;
		H5Fclose(hdf5_file);
		return false;
	}

	hid_t dataspace_id = H5Screate_simple(1, &size, NULL);

	/* Create a dataset attribute. */
	hid_t attribute_id = H5Acreate(loc, attr_name.c_str(), mem_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	if (attribute_id<0)
	{
		cerr << "HDF5_File_Writer::WriteAtrribute: Error, failed to create the attrbute" << endl;
		H5Sclose(dataspace_id);
		H5Oclose(loc);
		H5Fclose(hdf5_file);
		return false;
	}

	/* Write the attribute data. */
	if (H5Awrite(attribute_id, mem_type, value)<0)
	{
		cerr << "HDF5_File_Writer::WriteAtrribute: Error, failed to write the attrbute" << endl;
		H5Aclose(attribute_id);
		H5Sclose(dataspace_id);
		H5Oclose(loc);
		H5Fclose(hdf5_file);
		return false;
	}
	H5Aclose(attribute_id);
	H5Sclose(dataspace_id);
	H5Oclose(loc);
	H5Fclose(hdf5_file);
	return true;
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, float const* value, hsize_t size)
{
	return WriteAtrribute(locName,attr_name, value, size, H5T_NATIVE_FLOAT);
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, double const* value, hsize_t size)
{
	return WriteAtrribute(locName,attr_name, value, size, H5T_NATIVE_DOUBLE);
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, vector<float> values)
{
	float* val = new float[values.size()];
	for (size_t n=0;n<values.size();++n)
		val[n]=values.at(n);
	bool ok = HDF5_File_Writer::WriteAtrribute(locName, attr_name,val,values.size(),H5T_NATIVE_FLOAT);
	delete[] val; val=NULL;
	return ok;
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, vector<double> values)
{
	double* val = new double[values.size()];
	for (size_t n=0;n<values.size();++n)
		val[n]=values.at(n);
	bool ok =  HDF5_File_Writer::WriteAtrribute(locName, attr_name, val, values.size(), H5T_NATIVE_DOUBLE);
	delete[] val; val=NULL;
	return ok;
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, float value)
{
	return HDF5_File_Writer::WriteAtrribute(locName, attr_name,&value,1, H5T_NATIVE_FLOAT);
}

bool HDF5_File_Writer::WriteAtrribute(std::string locName, std::string attr_name, double value)
{
	return HDF5_File_Writer::WriteAtrribute(locName, attr_name,&value,1, H5T_NATIVE_DOUBLE);
}
