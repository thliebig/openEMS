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
	m_hdf5_file = -1;
	//suppress hdf5 error output
	//H5Eset_auto(NULL, NULL);

	if (!ReadAttribute<double>("/", "openEMS_HDF5_version", m_file_ver, true))
		m_file_ver = 0;

	if (m_file_ver>0.2)
	{
		if (!ReadAttribute<bool>("/", "openEMS_legacy_format", m_legacy_fmt, true))
			m_legacy_fmt = false; // set as default/fallback
	}
	else if (m_file_ver>0)
		m_legacy_fmt = true; // openEMS version 0.2 (and older?) always used the legacy format
}

hid_t HDF5_File_Reader::getFile()
{
	if (m_hdf5_file>=0)
		return m_hdf5_file;
	open();
	return m_hdf5_file;
}

bool HDF5_File_Reader::open()
{
	if (m_hdf5_file>=0)
		return true;
	m_hdf5_file = H5Fopen( m_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
	if (m_hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::open: opening the given file """ << m_filename << """ failed" << endl;
		return 0;
	}
	return (m_hdf5_file>=0);
}

void HDF5_File_Reader::close()
{
	if (m_hdf5_file<0)
		return;
	H5Fclose(m_hdf5_file);
	m_hdf5_file = -1;
}


HDF5_File_Reader::~HDF5_File_Reader()
{
	close();
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

bool HDF5_File_Reader::OpenGroup(hid_t &file, hid_t &group, string groupName)
{
	file = getFile();
	if (file==-1)
	{
		cerr << "HDF5_File_Reader::OpenGroup: opening the given file """ << m_filename << """ failed" << endl;
		return 0;
	}
	if (H5Lexists(file, groupName.c_str(), H5P_DEFAULT)<=0)
	{
		return 0;
	}

	group = H5Gopen(file, groupName.c_str(), H5P_DEFAULT );
	if (group<0)
	{
		cerr << "HDF5_File_Reader::OpenGroup: can't open group """ << groupName << """" << endl;
		return 0;
	}
	return true;
}

template <typename T>
bool HDF5_File_Reader::ReadAttribute(std::string grp_name, std::string attr_name, T &value, bool silent)
{
	hid_t hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}
	// H5Fclose(hdf5_file);
	// return false;

	if (H5Lexists(hdf5_file, grp_name.c_str(), H5P_DEFAULT)<=0)
	{
		if (!silent)
			cerr << "HDF5_File_Reader::ReadAttribute: group with the given name """ << grp_name << """ does not exist" << endl;
		return false;
	}

	if (H5Aexists_by_name(hdf5_file, grp_name.c_str(), attr_name.c_str(), H5P_DEFAULT)<=0)
	{
		if (!silent)
			cerr << "HDF5_File_Reader::ReadAttribute: attribue with the given name """ << attr_name << """ does not exist" << endl;
		return 0;
	}

	hid_t attr = H5Aopen_by_name(hdf5_file, grp_name.c_str(), attr_name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
	if (attr==-1)
	{
		if (!silent)
			cerr << "HDF5_File_Reader::ReadAttribute: Opening the given Attribute: """ << attr_name << """ failed" << endl;
		return false;
	}

	hid_t type = H5Aget_type(attr);
	if (type<0)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Dataset type error" << endl;
		H5Aclose(attr);
		return false;
	}

	size_t numVal = H5Aget_storage_size(attr)/H5Tget_size(type);
	if (numVal!=1)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Invalid attribute data size" << endl;
		H5Aclose(attr);
		return false;
	}

	if (!HDF5_File_Reader::CheckH5Type<T>(type))
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Unexpected attribute data type" << endl;
		H5Aclose(attr);
		return false;
	}
	if (H5Aread(attr, type, &value)<0)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Reading the given Attribute failed" << endl;
		H5Aclose(attr);
		return false;
	}

	H5Aclose(attr);
	return true;
}

bool HDF5_File_Reader::ReadAttribute(string grp_name, string attr_name, vector<float> &attr_values)
{
	vector<double> d_attr_values;
	if (ReadAttribute(grp_name, attr_name, d_attr_values)==false)
		return false;
	attr_values.resize(d_attr_values.size(),0);
	for (size_t n=0;n<d_attr_values.size();++n)
		attr_values.at(n)=d_attr_values.at(n);
	return true;
}

bool HDF5_File_Reader::ReadAttribute(string grp_name, string attr_name, vector<double> &attr_values)
{
	attr_values.clear();

	hid_t hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: opening the given file """ << m_filename << """ failed" << endl;
		return 0;
	}

	if (H5Lexists(hdf5_file, grp_name.c_str(), H5P_DEFAULT)<=0)
	{
		return false;
	}

	hid_t attr = H5Aopen_by_name(hdf5_file, grp_name.c_str(), attr_name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
	if (attr==-1)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Opening the given Attribute: """ << attr_name << """ failed" << endl;
		return false;
	}

	hid_t type = H5Aget_type(attr);
	if (type<0)
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Dataset type error" << endl;
		H5Aclose(attr);
		return false;
	}

	attr_values.clear();

	if (H5Tget_class(type)==H5T_FLOAT)
	{
		size_t numVal = H5Aget_storage_size(attr)/H5Tget_size(type);
		hid_t datatype=-1;
		void *value=NULL;
		float *f_value=NULL;
		double *d_value=NULL;
		if (H5Tget_size(type)==sizeof(float))
		{
			f_value = new float[numVal];
			value = f_value;
			datatype = H5T_NATIVE_FLOAT;
		}
		if (H5Tget_size(type)==sizeof(double))
		{
			d_value = new double[numVal];
			value = d_value;
			datatype = H5T_NATIVE_DOUBLE;
		}
		if (H5Aread(attr, datatype, value)<0)
		{
		  cerr << "HDF5_File_Reader::ReadAttribute: Reading the given Attribute failed" << endl;
		  H5Aclose(attr);
		  return false;
		}
		if (f_value)
			for (size_t n=0;n<numVal;++n)
				attr_values.push_back(f_value[n]);
		if (d_value)
			for (size_t n=0;n<numVal;++n)
				attr_values.push_back(d_value[n]);
		delete[] f_value;
		delete[] d_value;
	}
	else
	{
		cerr << "HDF5_File_Reader::ReadAttribute: Attribute type not supported" << endl;
		H5Aclose(attr);
		return false;
	}
	H5Aclose(attr);
	return true;
}

bool HDF5_File_Reader::ReadDataSet(string ds_name, hsize_t &nDim, hsize_t* &dims, float* &data)
{
	float* f_data;
	double* d_data;
	if (!ReadDataSet(ds_name, nDim, dims, f_data, d_data))
		return false;
	if (f_data != NULL)
	{
		data = f_data;
		return true;
	}

	hsize_t data_size = 1;
	for (unsigned int d=0;d<nDim;++d)
		data_size*=dims[d];
	data = new float[data_size];
	for (size_t n=0;n<data_size;++n)
		data[n]=d_data[n];
	delete[] d_data;
	d_data = NULL;
	return true;
}

bool HDF5_File_Reader::ReadDataSet(string ds_name, hsize_t &nDim, hsize_t* &dims, double* &data)
{
	float* f_data;
	double* d_data;
	if (!ReadDataSet(ds_name, nDim, dims, f_data, d_data))
		return false;
	if (d_data != NULL)
	{
		data = d_data;
		return true;
	}

	hsize_t data_size = 1;
	for (unsigned int d=0;d<nDim;++d)
		data_size*=dims[d];
	data = new double[data_size];
	for (size_t n=0;n<data_size;++n)
		data[n]=f_data[n];
	delete[] f_data;
	f_data = NULL;
	return true;
}

template <typename T>
bool HDF5_File_Reader::ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, T* &data)
{
	return false;
}

template <typename T>
hid_t HDF5_File_Reader::GetH5Type()
{
	if (typeid(T) == typeid(float))
		return H5T_NATIVE_FLOAT;
	else if (typeid(T) == typeid(double))
		return H5T_NATIVE_DOUBLE;
	else if (typeid(T) == typeid(int))
		return H5T_NATIVE_INT;
	else if (typeid(T) == typeid(unsigned int))
		return H5T_NATIVE_UINT;
	else if (typeid(T) == typeid(bool))
		return H5T_NATIVE_HBOOL;
	else if (typeid(T) == typeid(complex<float>))
	{
		typedef struct {
		float re;   /*real part */
		float im;   /*imaginary part */
		} complex_t;

		// create numpy compatible complex128
		complex_t tmp;  /*used only to compute offsets */
		hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof tmp);
		H5Tinsert (complex_id, "r", HOFFSET(complex_t,re), H5T_NATIVE_FLOAT);
		H5Tinsert (complex_id, "i", HOFFSET(complex_t,im), H5T_NATIVE_FLOAT);
		return complex_id;
	}
	else if (typeid(T) == typeid(complex<double>))
	{
		typedef struct {
		double re;   /*real part */
		double im;   /*imaginary part */
		} complex_t;

		// create numpy compatible complex128
		complex_t tmp;  /*used only to compute offsets */
		hid_t complex_id = H5Tcreate (H5T_COMPOUND, sizeof tmp);
		H5Tinsert (complex_id, "r", HOFFSET(complex_t,re), H5T_NATIVE_DOUBLE);
		H5Tinsert (complex_id, "i", HOFFSET(complex_t,im), H5T_NATIVE_DOUBLE);
		return complex_id;
	}
	return H5I_INVALID_HID; // H5I_INVALID_HID should be -1
}

// Explicit template instantiation
template hid_t HDF5_File_Reader::GetH5Type<int>();
template hid_t HDF5_File_Reader::GetH5Type<float>();
template hid_t HDF5_File_Reader::GetH5Type<double>();
template hid_t HDF5_File_Reader::GetH5Type<unsigned int>();
template hid_t HDF5_File_Reader::GetH5Type<bool>();
template hid_t HDF5_File_Reader::GetH5Type<complex<float>>();
template hid_t HDF5_File_Reader::GetH5Type<complex<double>>();

template <typename T>
bool HDF5_File_Reader::CheckH5Type(hid_t type)
{
	// cerr << "CheckH5Type Comp?: " << (H5Tget_class(type)==H5T_COMPOUND) << " size? " << H5Tget_size(type) << endl;
	if (H5Tequal(type, H5T_NATIVE_FLOAT) && (typeid(T) == typeid(float)))
		return true;
	else if (H5Tequal(type, H5T_NATIVE_DOUBLE) && (typeid(T) == typeid(double)))
		return true;
	else if (H5Tequal(type, H5T_NATIVE_INT) && (typeid(T) == typeid(int)))
		return true;
	else if (H5Tequal(type, H5T_NATIVE_UINT) && (typeid(T) == typeid(unsigned int)))
		return true;
	else if (H5Tequal(type, H5T_NATIVE_HBOOL) && (typeid(T) == typeid(bool)))
		return true;
	else if ((typeid(T) == typeid(complex<float>)) && (H5Tget_class(type)==H5T_COMPOUND) && (H5Tget_size(type)==2*sizeof(float)))
		return true;
	else if ((typeid(T) ==typeid(complex<double>)) && (H5Tget_class(type)==H5T_COMPOUND) && (H5Tget_size(type)==2*sizeof(double)))
		return true;
	return false;
}

bool HDF5_File_Reader::GetDataSetInfo(std::string ds_name, hid_t &type, hsize_t &nDim, hsize_t* &dims)
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hid_t dataset = H5Dopen(hdf5_file, ds_name.c_str(), H5P_DEFAULT );
	if (dataset<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: dataset not found" << endl;
		return false;
	}
	type = H5Dget_type(dataset);

	hid_t space = H5Dget_space(dataset);
	nDim = H5Sget_simple_extent_ndims(space);
	dims = new hsize_t[nDim];
	H5Sget_simple_extent_dims(space, dims, NULL );

	H5Dclose(dataset);
	return true;
}

bool HDF5_File_Reader::CheckDataSetExists(std::string ds_name)
{
	hid_t hdf5_file = getFile();
	if (hdf5_file<0)
		return false;

	bool ok = H5Lexists(hdf5_file, ds_name.c_str(), H5P_DEFAULT);
	return ok;
}

bool HDF5_File_Reader::ReadDataSet(std::string ds_name, hid_t type, void* buffer)
{
	hid_t hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hid_t dataset = H5Dopen(hdf5_file, ds_name.c_str(), H5P_DEFAULT );
	if (dataset<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: dataset not found" << endl;
		return false;
	}

	if (H5Dread(dataset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,buffer)<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: error reading data" << endl;
		H5Dclose(dataset);
		return false;
	}
	H5Dclose(dataset);
	return true;
}


bool HDF5_File_Reader::ReadDataSet(std::string ds_name, hsize_t &nDim, hsize_t* &dims, float* &f_data, double* &d_data)
{
	hid_t type;
	if (GetDataSetInfo(ds_name, type, nDim, dims)==false)
		return false;

	hid_t hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}

	hid_t dataset = H5Dopen(hdf5_file, ds_name.c_str(), H5P_DEFAULT );
	if (dataset<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: dataset not found" << endl;
		return false;
	}

	hsize_t data_size = 1;
	for (unsigned int d=0;d<nDim;++d)
		data_size*=dims[d];

	void *value=NULL;
	if (H5Tequal(type,H5T_NATIVE_FLOAT))
	{
		f_data = new float[data_size];
		d_data = NULL;
		value = f_data;
	}
	else if (H5Tequal(type,H5T_NATIVE_DOUBLE))
	{
		f_data = NULL;
		d_data = new double[data_size];
		value = d_data;
	}
	else
	{
		H5Dclose(dataset);
		f_data = NULL;
		d_data = NULL;
		return false;
	}

	if (H5Dread(dataset,type,H5S_ALL,H5S_ALL,H5P_DEFAULT,value)<0)
	{
		cerr << "HDF5_File_Reader::ReadDataSet: error reading data" << endl;
		H5Dclose(dataset);
		delete[] f_data;
		f_data = NULL;
		delete[] d_data;
		d_data = NULL;
		return false;
	}
	H5Dclose(dataset);
	return true;
}

bool HDF5_File_Reader::ReadMesh(float** lines, unsigned int* numLines, int &meshType, std::string s_mesh_grp)
{
	double* d_lines[3]={NULL,NULL,NULL};
	if (!ReadMesh(d_lines, numLines, meshType, s_mesh_grp))
		return false;
	for (int n=0;n<3;++n)
	{
		lines[n] = new float[numLines[n]];
		for (unsigned int i=0; i<numLines[n]; ++i)
			lines[n][i]=d_lines[n][i];
	}
	return true;
}

bool HDF5_File_Reader::ReadMesh(double** lines, unsigned int* numLines, int &meshType, std::string s_mesh_grp)
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file = getFile();
	if (hdf5_file==-1)
	{
		cerr << "HDF5_File_Reader::ReadMesh: opening the given file """ << m_filename << """ failed" << endl;
		return false;
	}
	std::string sx = s_mesh_grp + "/x";
	std::string sr = s_mesh_grp + "/rho";
	std::string sy = s_mesh_grp + "/y";
	std::string sa = s_mesh_grp + "/alpha";
	std::string sz = s_mesh_grp + "/z";
	vector<string> Mesh_Names;
	if (H5Lexists(hdf5_file, sx.c_str(), H5P_DEFAULT) && H5Lexists(hdf5_file, sy.c_str(), H5P_DEFAULT) && H5Lexists(hdf5_file, sz.c_str(), H5P_DEFAULT))
	{
		meshType = 0;
		Mesh_Names.push_back(sx);
		Mesh_Names.push_back(sy);
		Mesh_Names.push_back(sz);
	}
	else if (H5Lexists(hdf5_file, sr.c_str(), H5P_DEFAULT) && H5Lexists(hdf5_file, sa.c_str(), H5P_DEFAULT) && H5Lexists(hdf5_file, sz.c_str(), H5P_DEFAULT))
	{
		meshType = 1;
		Mesh_Names.push_back(sr);
		Mesh_Names.push_back(sa);
		Mesh_Names.push_back(sz);
	}
	else
	{
		cerr << "HDF5_File_Reader::ReadMesh: no falid mesh information found" << endl;
		return false;
	}

	for (int n=0;n<3;++n)
	{
		hsize_t nDim;
		hsize_t* dims=NULL;
		double* data=NULL;
		ReadDataSet(Mesh_Names.at(n), nDim, dims, data);
		if (nDim!=1)
		{
			cerr << "HDF5_File_Reader::ReadMesh: mesh dimension error" << endl;
			delete[] dims;
			delete[] data;
			return false;
		}
		numLines[n]=dims[0];
		delete[] dims;
		lines[n]=data;
	}
	return true;
}

unsigned int HDF5_File_Reader::GetNumTimeSteps()
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file;
	hid_t TD_grp;
	if (OpenGroup(hdf5_file, TD_grp, "/FieldData/TD")==false)
		return false;

	hsize_t numObj;
	if (H5Gget_num_objs(TD_grp,&numObj)<0)
	{
		cerr << "HDF5_File_Reader::GetNumTimeSteps: can't read number of datasets" << endl;
		H5Gclose(TD_grp);
		return 0;
	}
	H5Gclose(TD_grp);
	return numObj;
}

bool HDF5_File_Reader::ReadTimeSteps(vector<unsigned int> &timestep, vector<string> &names)
{
	if (IsValid()==false)
		return false;

	hid_t hdf5_file;
	hid_t TD_grp;
	if (OpenGroup(hdf5_file, TD_grp, "/FieldData/TD")==false)
		return false;

	hsize_t numObj;
	if (H5Gget_num_objs(TD_grp,&numObj)<0)
	{
		cerr << "HDF5_File_Reader::ReadTimeSteps: can't read number of datasets" << endl;
		H5Gclose(TD_grp);
		return false;
	}

	std::string name;
	timestep.clear();
	timestep.resize(numObj,0);
	names.clear();
	names.resize(numObj);
	for (hsize_t n=0;n<numObj;++n)
	{
		if (!GetDataSetNameByIndex(TD_grp, n, name))
		{
			cerr << "HDF5_File_Reader::ReadTimeSteps: invalid timestep found!" << endl;
			H5Gclose(TD_grp);
			return false;
		}

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
			return false;
		}
	}
	H5Gclose(TD_grp);
	return true;
}

bool HDF5_File_Reader::GetDataSetNameByIndex(std::string grp_name, unsigned int idx, std::string &name)
{
	if (IsValid()==false)
		return false;
	hid_t grp;
	hid_t file;
	if (!OpenGroup(file, grp, grp_name))
		return false;

	bool ok = GetDataSetNameByIndex(grp, idx, name);
	H5Gclose(grp);
	return ok;
}

bool HDF5_File_Reader::GetDataSetNameByIndex(hid_t &group, unsigned int idx, std::string &name)
{
	hsize_t numObj;
	if (H5Gget_num_objs(group,&numObj)<0)
	{
		cerr << "HDF5_File_Reader::GetDataSetNameByIndex: can't read number of datasets" << endl;
		return false;
	}

	if (idx>=numObj)
		return false;

	ssize_t len = H5Gget_objname_by_idx(group, idx, NULL, 0 );
	if (len<=0)
		return false;

	char c_name[len+1];
	len = H5Gget_objname_by_idx(group, idx, c_name, len+1);
	name = c_name;
	return true;
}


bool HDF5_File_Reader::GetTDVectorData(size_t idx, float &time, ArrayLib::ArrayNIJK<float> &data)
{
	if (IsValid()==false)
		return false;

	std::string ds_name;
	if (!GetDataSetNameByIndex("/FieldData/TD", idx, ds_name))
		return false;

	ds_name = "/FieldData/TD/" + ds_name;
	if (!ReadAttribute<float>(ds_name, "time", time, true))
	{
		cerr << "HDF5_File_Reader::GetTDVectorData: can't read time attribute!" << endl;
		return false;
	}
	return ReadVectorDataSet<float>(ds_name, data);
}

unsigned int HDF5_File_Reader::GetNumFrequencies()
{
	vector<float> frequencies;
	if (ReadFrequencies(frequencies)==false)
		return 0;
	return frequencies.size();
}

bool HDF5_File_Reader::ReadFrequencies(vector<float> &frequencies)
{
	if (IsValid()==false)
		return false;

	return ReadAttribute("/FieldData/FD","frequency",frequencies);
}

bool HDF5_File_Reader::ReadFrequencies(vector<double> &frequencies)
{
	if (IsValid()==false)
		return false;

	return ReadAttribute("/FieldData/FD","frequency",frequencies);
}

template <typename T>
bool HDF5_File_Reader::ReadScalarDataSet(std::string ds_name, ArrayLib::ArrayIJK<T> &data, bool type_check)
{
	if (IsValid()==false)
		return false;

	hsize_t nDim;
	hsize_t* dims=NULL;
	hid_t type;

	if (!GetDataSetInfo(ds_name, type, nDim, dims))
	{
		return false;
	}

	if ((type_check) && (!HDF5_File_Reader::CheckH5Type<T>(type)))
	{
		cerr << "HDF5_File_Reader::ReadVectorDataSet: Unexpected dataset data type" << endl;
		return false;
	}

	// make sure all dimensions are as expected!
	if ((nDim!=3))
	{
		delete[] dims;
		return false;
	}
	unsigned int ui_dim[3];
	for (int n=0;n<3;++n)
		ui_dim[n] = dims[n];

	data.Init(ds_name, ui_dim);
	bool ok = ReadDataSet(ds_name, type, (void*)data.data());
	delete[] dims;
	return ok;

	return false;
}

// Explicit template instantiation
template bool HDF5_File_Reader::ReadScalarDataSet<float>(std::string ds_name, ArrayLib::ArrayIJK<float> &data, bool type_check);

bool HDF5_File_Reader::GetFDVectorData(size_t idx, ArrayLib::ArrayNIJK<std::complex<float>> &data)
{
	stringstream ds_name;
	ds_name << "/FieldData/FD/f" << idx;
	return ReadVectorDataSet<std::complex<float>>(ds_name.str(), data);
}

template <typename T>
bool HDF5_File_Reader::ReadVectorDataSet(std::string ds_name, ArrayLib::ArrayNIJK<T> &data, bool type_check)
{
	if (IsValid()==false)
		return false;

	hsize_t nDim;
	hsize_t* dims=NULL;
	hid_t type;

	if (!GetDataSetInfo(ds_name, type, nDim, dims))
	{
		return false;
	}

	if ((type_check) && (!HDF5_File_Reader::CheckH5Type<T>(type)))
	{
		cerr << "HDF5_File_Reader::ReadVectorDataSet: Unexpected dataset data type" << endl;
		return false;
	}

	// make sure all dimensions are as expected!
	if ((nDim!=4) || (dims[0]!=3))
	{
		delete dims;
		return false;
	}

	unsigned int ui_dim[3];
	for (int n=0;n<3;++n)
		ui_dim[n] = dims[n+1];


	data.Init(ds_name, ui_dim);
	bool ok = ReadDataSet(ds_name, type, (void*)data.data());
	delete[] dims;
	return ok;

	return false;
}

// Explicit template instantiation
template bool HDF5_File_Reader::ReadVectorDataSet<float>(std::string ds_name, ArrayLib::ArrayNIJK<float> &data, bool type_check);


bool HDF5_File_Reader::CalcFDVectorData(vector<float> &frequencies, std::vector<ArrayLib::ArrayNIJK<std::complex<float>>*> &FD_data)
{
	for (unsigned int n=0;n<FD_data.size();++n)
		delete FD_data.at(n);
	FD_data.clear();

	unsigned int numTS = GetNumTimeSteps();
	if (numTS<=0)
	{
		cerr << "HDF5_File_Reader::CalcFDVectorData: error, no TD data found..." << endl;
		return false;
	}

	float time;
	//read first TD data
	ArrayLib::ArrayNIJK<float> field;
	if (!GetTDVectorData(0,time,field))
	{
		cerr << "HDF5_File_Reader::CalcFDVectorData: error, no TD data found..." << endl;
		return false;
	}

	//init
	unsigned int datasize[3] = {field.extent(1), field.extent(2), field.extent(3)};
	FD_data.resize(frequencies.size(), NULL);
	for (size_t fn=0;fn<frequencies.size();++fn)
	{
		ArrayLib::ArrayNIJK<std::complex<float>>* field_arr = new ArrayLib::ArrayNIJK<std::complex<float>>("FD", datasize);
		FD_data.at(fn) = field_arr;
	}

	size_t ts=0;
	unsigned int pos[3];
	complex<float> PI_2_I(0.0,-2.0*M_PI);
	complex<float> exp_jwt_2_dt;
	float time_diff=0;
	float time_old =0;
	size_t N = field.size();
	std::complex<float>* field_fd;
	float* td_data;
	// field.Reset();
	while (field.valid())
	{
		if ((ts>1) && abs(time_diff - (time - time_old))>1e15)
		{
			cerr << "HDF5_File_Reader::CalcFDVectorData: time interval error..." << endl;
			for (size_t fn=0;fn<frequencies.size();++fn)
				delete FD_data.at(fn);
			FD_data.clear();
			return false;
		}
		time_diff = time - time_old;
		td_data = field.data();
		for (size_t fn=0;fn<frequencies.size();++fn)
		{
			exp_jwt_2_dt = exp( (complex<float>)(PI_2_I * frequencies.at(fn) * time) );
			field_fd = FD_data.at(fn)->data();
			for (size_t n=0;n<N;++n)
				field_fd[n] += td_data[n] * exp_jwt_2_dt;
		}
		++ts;
		field.Reset();
		time_old = time;
		GetTDVectorData(ts,time,field);
	}

	// finalize data
	time_diff*=2;
	for (size_t fn=0;fn<frequencies.size();++fn)
	{
		field_fd = FD_data.at(fn)->data();
		for (size_t n=0;n<N;++n)
			field_fd[n] *= time_diff;
	}
	return true;
}
