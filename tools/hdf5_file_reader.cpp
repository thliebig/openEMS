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

#include "hdf5_file_reader.h"
#include "../tools/array_ops.h"
#include <hdf5.h>

#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;

HDF5_File_Reader::HDF5_File_Reader(string filename)
{
	m_filename = filename;
	//suppress hdf5 error output
	//H5Eset_auto(NULL, NULL);
}

HDF5_File_Reader::~HDF5_File_Reader()
{
}

bool HDF5_File_Reader::IsValid()
{
	htri_t val = H5Fis_hdf5(m_filename.c_str());
	if (val<0)
	{
		cerr << "HDF5_File_Reader::IsValid: the given file """ << m_filename << """ is not accessible..." << endl;
		return false;
	}
	if (val==0)
	{
		cerr << "HDF5_File_Reader::IsValid: the given file """ << m_filename << """ is invalid..." << endl;
		return false;
	}
	if (val==0)
	cerr << "HDF5_File_Reader::IsValid: the given file """ << m_filename << """ is valid..." << endl;
	return true;
}

bool HDF5_File_Reader::ReadDataSet(string ds_name, hsize_t &nDim, hsize_t* &dims, float* &data)
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hid_t dataset = H5Dopen(hdf5_file, ds_name.c_str() );
	if (dataset<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: dataset not found" << endl;
		H5Fclose(hdf5_file);
		return false;
	}
	hid_t type = H5Dget_type(dataset);
	if (type<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: dataset type error" << endl;
		H5Dclose(dataset);
		H5Fclose(hdf5_file);
		return false;
	}
	if (H5Tget_class(type)!=H5T_FLOAT)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: dataset type not a float" << endl;
		H5Dclose(dataset);
		H5Fclose(hdf5_file);
		return false;
	}

	hid_t space = H5Dget_space(dataset);
	nDim = H5Sget_simple_extent_ndims(space);
	dims = new hsize_t[nDim];
	H5Sget_simple_extent_dims(space, dims, NULL );
	hsize_t data_size = 1;
	for (unsigned int d=0;d<nDim;++d)
		data_size*=dims[d];
	data = new float[data_size];
	if (H5Dread(dataset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,data)<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: error reading data" << endl;
		H5Dclose(dataset);
		H5Fclose(hdf5_file);
		delete[] data;
		data=NULL;
		return false;
	}
	H5Dclose(dataset);
	H5Fclose(hdf5_file);
	return true;
}

bool HDF5_File_Reader::ReadMesh(float** lines, unsigned int* numLines, int &meshType)
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadMesh: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	vector<string> Mesh_Names;
	if (H5Lexists(hdf5_file, "/Mesh/x", H5P_DEFAULT) && H5Lexists(hdf5_file, "/Mesh/y", H5P_DEFAULT) && H5Lexists(hdf5_file, "/Mesh/z", H5P_DEFAULT))
	{
		meshType = 0;
		Mesh_Names.push_back("/Mesh/x");
		Mesh_Names.push_back("/Mesh/y");
		Mesh_Names.push_back("/Mesh/z");
	}
	else if (H5Lexists(hdf5_file, "/Mesh/rho", H5P_DEFAULT) && H5Lexists(hdf5_file, "/Mesh/alpha", H5P_DEFAULT) && H5Lexists(hdf5_file, "/Mesh/z", H5P_DEFAULT))
	{
		meshType = 1;
		Mesh_Names.push_back("/Mesh/rho");
		Mesh_Names.push_back("/Mesh/alpha");
		Mesh_Names.push_back("/Mesh/z");
	}
	else
	{
		cerr << "HDF5_File_Reader::ReadMesh: no falid mesh information found" << endl;
		H5Fclose(hdf5_file);
		return false;
	}
	for (int n=0;n<3;++n)
	{
		hsize_t nDim;
		hsize_t* dims=NULL;
		float* data=NULL;
		ReadDataSet(Mesh_Names.at(n), nDim, dims, data);
		if (nDim!=1)
		{
			cerr << "HDF5_File_Reader::ReadMesh: mesh dimension error" << endl;
			delete[] dims;
			delete[] data;
			H5Fclose(hdf5_file);
			return false;
		}
		numLines[n]=dims[0];
		delete[] dims;
		lines[n]=data;
	}
	H5Fclose(hdf5_file);
	return true;
}

unsigned int HDF5_File_Reader::GetNumTimeSteps()
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::GetNumTimeSteps: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	if (H5Lexists(hdf5_file, "/FieldData/TD", H5P_DEFAULT)<=0)
	{
		H5Fclose(hdf5_file);
		return 0;
	}

	hid_t TD_grp = H5Gopen(hdf5_file, "/FieldData/TD" );
	if (TD_grp<0)
	{
		cerr << "HDF5_File_Reader::GetNumTimeSteps: can't open group ""/FieldData/TD""" << endl;
		H5Fclose(hdf5_file);
		return 0;
	}
	hsize_t numObj;
	if (H5Gget_num_objs(TD_grp,&numObj)<0)
	{
		cerr << "HDF5_File_Reader::GetNumTimeSteps: can't read number of datasets" << endl;
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return 0;
	}
	H5Gclose(TD_grp);
	H5Fclose(hdf5_file);
	return numObj;
}

bool HDF5_File_Reader::ReadTimeSteps(vector<unsigned int> &timestep, vector<string> &names)
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadTimeSteps: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	if (H5Lexists(hdf5_file, "/FieldData/TD", H5P_DEFAULT)<0)
	{
		cerr << "HDF5_File_Reader::ReadTimeSteps: can't open ""/FieldData/TD""" << endl;
		H5Fclose(hdf5_file);
		return false;
	}
	hid_t TD_grp = H5Gopen(hdf5_file, "/FieldData/TD" );
	if (TD_grp<0)
	{
		cerr << "HDF5_File_Reader::ReadTimeSteps: can't open ""/FieldData/TD""" << endl;
		H5Fclose(hdf5_file);
		return false;
	}
	hsize_t numObj;
	if (H5Gget_num_objs(TD_grp,&numObj)<0)
	{
		cerr << "HDF5_File_Reader::ReadTimeSteps: can't read number of datasets" << endl;
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return false;
	}
	char name[100];
	timestep.clear();
	timestep.resize(numObj,0);
	names.clear();
	names.resize(numObj);
	for (hsize_t n=0;n<numObj;++n)
	{
		if (H5Gget_objtype_by_idx(TD_grp, n)  != H5G_DATASET)
		{
			cerr << "HDF5_File_Reader::ReadTimeSteps: invalid timestep found!" << endl;
			H5Gclose(TD_grp);
			H5Fclose(hdf5_file);
			return false;
		}
		ssize_t name_size = H5Gget_objname_by_idx(TD_grp, n, name, 100 );
		istringstream is(name);
		unsigned int num;
		if (is >> num)
		{
			timestep.at(n)=num;
			names.at(n)=name;
		}
		else
		{
			cerr << "HDF5_File_Reader::ReadTimeSteps: invalid timestep format found!" << endl;
			H5Gclose(TD_grp);
			H5Fclose(hdf5_file);
			return false;
		}
	}
	H5Gclose(TD_grp);
	H5Fclose(hdf5_file);
	return true;
}

float**** HDF5_File_Reader::GetTDVectorData(size_t idx, float &time, unsigned int data_size[])
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	time = 0;
	if (H5Lexists(hdf5_file, "/FieldData/TD", H5P_DEFAULT)<=0)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: can't open ""/FieldData/TD""" << endl;
		H5Fclose(hdf5_file);
		return NULL;
	}
	hid_t TD_grp = H5Gopen(hdf5_file, "/FieldData/TD" );
	if (TD_grp<0)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: can't open ""/FieldData/TD""" << endl;
		H5Fclose(hdf5_file);
		return NULL;
	}

	hsize_t numObj;
	if (H5Gget_num_objs(TD_grp,&numObj)<0)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: can't read number of datasets" << endl;
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}
	if (idx>=numObj)
	{
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}

	if (H5Gget_objtype_by_idx(TD_grp, idx)  != H5G_DATASET)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: invalid timestep found!" << endl;
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}

	char name[100];
	ssize_t name_size = H5Gget_objname_by_idx(TD_grp, idx, name, 100 );
	string ds_name = "/FieldData/TD/" + string(name);

	hid_t attr = H5Aopen_by_name(hdf5_file, ds_name.c_str(), "time", H5P_DEFAULT, H5P_DEFAULT);
	if (attr<0)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: time attribute not found!" << endl;
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}
	if (H5Aread(attr, H5T_NATIVE_FLOAT, &time)<0)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: can't read time attribute!" << endl;
		H5Aclose(attr);
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}

	hsize_t nDim;
	hsize_t* dims=NULL;
	float* data=NULL;
	ReadDataSet(ds_name, nDim, dims, data);
	if (nDim!=4)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: data dimension invalid" << endl;
		delete[] dims;
		H5Aclose(attr);
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}
	if (dims[0]!=3)
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: vector data dimension invalid" << endl;
		delete[] dims;
		H5Aclose(attr);
		H5Gclose(TD_grp);
		H5Fclose(hdf5_file);
		return NULL;
	}
	data_size[0]=dims[3];
	data_size[1]=dims[2];
	data_size[2]=dims[1];
	delete[] dims;
	data_size[3]=3;
	size_t pos = 0;
	float**** field = Create_N_3DArray<float>(data_size);
	for (unsigned int d=0;d<3;++d)
		for (unsigned int k=0;k<data_size[2];++k)
			for (unsigned int j=0;j<data_size[1];++j)
				for (unsigned int i=0;i<data_size[0];++i)
				{
					field[d][i][j][k]=data[pos++];
				}
	delete[] data;
	H5Aclose(attr);
	H5Gclose(TD_grp);
	H5Fclose(hdf5_file);
	return field;
}

unsigned int HDF5_File_Reader::GetNumFrequencies()
{
	cerr << "HDF5_File_Reader::GetNumFrequencies(): not implemented yet!" << endl;
	return 0;
}

bool HDF5_File_Reader::ReadFrequencies(vector<float> &frequencies)
{
	cerr << "HDF5_File_Reader::ReadFrequencies(): not implemented yet!" << endl;
	return false;
}

complex<float>**** HDF5_File_Reader::GetFDVectorData(size_t idx, float &frequency, unsigned int data_size[])
{
	cerr << "HDF5_File_Reader::GetFDVectorData(): not implemented yet!" << endl;
	return NULL;
}

bool HDF5_File_Reader::CalcFDVectorData(vector<float> &frequencies, vector<complex<float>****> &FD_data, unsigned int data_size[4])
{
	FD_data.clear();

	if (GetNumTimeSteps()<=0)
	{
		cerr << "HDF5_File_Reader::CalcFDVectorData: error, no TD data found..." << endl;
		return false;
	}

	float time;
	//read first TD data
	float**** field = this->GetTDVectorData(0,time,data_size);
	if (field==NULL)
	{
		cerr << "HDF5_File_Reader::CalcFDVectorData: error, no TD data found..." << endl;
		return false;
	}

	//init
	FD_data.resize(frequencies.size(), NULL);
	for (size_t fn=0;fn<frequencies.size();++fn)
		FD_data.at(fn) = Create_N_3DArray<complex<float> >(data_size);

	size_t ts=0;
	unsigned int pos[3];
	complex<float> PI_2_I(0.0,-2.0*M_PI);
	complex<float> exp_jwt_2_dt;
	float time_diff=0;
	float time_old =0;
	complex<float>**** field_fd = NULL;
	while (field)
	{
		if ((ts>1) && abs(time_diff - (time - time_old))>1e15)
		{
			cerr << "HDF5_File_Reader::CalcFDVectorData: time interval error..." << endl;
			for (size_t fn=0;fn<frequencies.size();++fn)
				Delete_N_3DArray(FD_data.at(fn),data_size);
			FD_data.clear();
			return false;
		}
		time_diff = time - time_old;
		for (size_t fn=0;fn<frequencies.size();++fn)
		{
			exp_jwt_2_dt = exp( (complex<float>)(PI_2_I * frequencies.at(fn) * time) );
			field_fd = FD_data.at(fn);
			for (pos[0]=0; pos[0]<data_size[0]; ++pos[0])
			{
				for (pos[1]=0; pos[1]<data_size[1]; ++pos[1])
				{
					for (pos[2]=0; pos[2]<data_size[2]; ++pos[2])
					{
						field_fd[0][pos[0]][pos[1]][pos[2]] += field[0][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
						field_fd[1][pos[0]][pos[1]][pos[2]] += field[1][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
						field_fd[2][pos[0]][pos[1]][pos[2]] += field[2][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
					}
				}
			}
		}
		++ts;
		Delete_N_3DArray(field,data_size);
		time_old = time;
		field = this->GetTDVectorData(ts,time,data_size);
	}

	// finalize data
	time_diff*=2;
	for (size_t fn=0;fn<frequencies.size();++fn)
	{
		field_fd = FD_data.at(fn);
		for (pos[0]=0; pos[0]<data_size[0]; ++pos[0])
		{
			for (pos[1]=0; pos[1]<data_size[1]; ++pos[1])
			{
				for (pos[2]=0; pos[2]<data_size[2]; ++pos[2])
				{
					field_fd[0][pos[0]][pos[1]][pos[2]] *= time_diff;
					field_fd[1][pos[0]][pos[1]][pos[2]] *= time_diff;
					field_fd[2][pos[0]][pos[1]][pos[2]] *= time_diff;
				}
			}
		}
	}
	return true;
}
