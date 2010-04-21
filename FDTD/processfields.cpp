/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "processfields.h"

#include <iomanip>
#include "H5Cpp.h"

ProcessFields::ProcessFields(Operator* op, Engine* eng) : Processing(op, eng)
{
	m_DumpMode = NO_INTERPOLATION;
	m_DumpType = E_FIELD_DUMP;
	// vtk-file is default
	m_fileType = VTK_FILETYPE;
	SetSubSampling(1);
	SetPrecision(6);

	for (int n=0;n<3;++n)
	{
		numLines[n]=0;
		discDLines[n]=NULL;
		discLines[n]=NULL;
	}
}

ProcessFields::~ProcessFields()
{
	for (int n=0;n<3;++n)
	{
		delete[] discDLines[n];
		discDLines[n]=NULL;
		delete[] discLines[n];
		discLines[n]=NULL;
	}
}

void ProcessFields::InitProcess()
{
	if (Enabled==false) return;
	string names[] = {"x","y","z"};
	if (m_fileType==HDF5_FILETYPE)
	{
		unsigned int* NrLines;
		double** Lines;

		if (m_DumpMode==CELL_INTERPOLATE)
		{
			NrLines = numDLines;
			Lines = discDLines;
		}
		else if (m_DumpMode==NO_INTERPOLATION)
		{
			NrLines = numLines;
			Lines = discLines;
		}
		else
			return;

		m_fileName+= ".h5";
		H5::H5File* file = new H5::H5File( m_fileName , H5F_ACC_TRUNC );

		H5::Group* group = new H5::Group( file->createGroup( "/Mesh" ));
		for (int n=0;n<3;++n)
		{
			hsize_t dimsf[1];              // dataset dimensions
			dimsf[0] = NrLines[n];
			H5::DataSpace dataspace( 1, dimsf );
			H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );
			H5::DataSet dataset = group->createDataSet( names[n].c_str(), datatype, dataspace );
			//convert to float...
			float* array = new float[NrLines[n]];
			for (unsigned int i=0;i<NrLines[n];++i)
				array[i] = Lines[n][i];
			//write to dataset
			dataset.write( array, H5::PredType::NATIVE_FLOAT );
		}
		delete group;

		group = new H5::Group( file->createGroup( "/FieldData" ));
		delete group;
		delete file;
	}
}

string ProcessFields::GetFieldNameByType(DumpType type)
{
	switch (type)
	{
		case E_FIELD_DUMP:
			return "E-Field";
		case H_FIELD_DUMP:
			return "H-Field";
	}
	return "unknown field";
}

void ProcessFields::DefineStartStopCoord(double* dstart, double* dstop)
{
	vector<double> lines;
	if (m_DumpMode==NO_INTERPOLATION)
	{
		if (Op->SnapToMesh(dstart,start)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (Op->SnapToMesh(dstop,stop)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;
		//create dual mesh
		for (int n=0;n<3;++n)
		{
	//		cerr << "start " << start[n] << "stop " << stop[n];
			if (start[n]>stop[n])
			{
				unsigned int help = start[n];
				start[n]=stop[n];
				stop[n]=help;
			}
			lines.clear();
			for (unsigned int i=start[n];i<=stop[n];i+=subSample[n])
			{
				lines.push_back(Op->GetDiscLine(n,i));//Op->discLines[n][i]);
			}
			numLines[n] = lines.size();
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			for (unsigned int i=0;i<numLines[n];++i)
				discLines[n][i] = lines.at(i);
		}
	}
	else if (m_DumpMode==CELL_INTERPOLATE)
	{
		if (Op->SnapToMesh(dstart,start,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (Op->SnapToMesh(dstop,stop,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;

		//create dual mesh
		for (int n=0;n<3;++n)
		{
	//		cerr << "start " << start[n] << "stop " << stop[n];
			if (start[n]>stop[n])
			{
				unsigned int help = start[n];
				start[n]=stop[n];
				stop[n]=help;
			}
			++stop[n];
			lines.clear();
			for (unsigned int i=start[n];i<stop[n];i+=subSample[n])
			{
				lines.push_back(Op->GetDiscLine(n,i,true));//0.5*(Op->discLines[n][i+1] +  Op->discLines[n][i]));
			}		
			numDLines[n] = lines.size();
			delete[] discDLines[n];
			discDLines[n] = new double[numDLines[n]];
			for (unsigned int i=0;i<numDLines[n];++i)
				discDLines[n][i] = lines.at(i);
		}
	}
}

double ProcessFields::CalcTotalEnergy() const
{
	if (!Eng)
		return 0;

	double energy=0;
	unsigned int pos[3];
	for (pos[0]=0;pos[0]<Op->GetNumberOfLines(0);++pos[0])
	{
		for (pos[1]=0;pos[1]<Op->GetNumberOfLines(1);++pos[1])
		{
			for (pos[2]=0;pos[2]<Op->GetNumberOfLines(2);++pos[2])
			{
				energy+=fabs(Eng->GetVolt(0,pos[0],pos[1],pos[2]) * Eng->GetCurr(1,pos[0],pos[1],pos[2]));
				energy+=fabs(Eng->GetVolt(0,pos[0],pos[1],pos[2]) * Eng->GetCurr(2,pos[0],pos[1],pos[2]));
				energy+=fabs(Eng->GetVolt(1,pos[0],pos[1],pos[2]) * Eng->GetCurr(0,pos[0],pos[1],pos[2]));
				energy+=fabs(Eng->GetVolt(1,pos[0],pos[1],pos[2]) * Eng->GetCurr(2,pos[0],pos[1],pos[2]));
				energy+=fabs(Eng->GetVolt(2,pos[0],pos[1],pos[2]) * Eng->GetCurr(0,pos[0],pos[1],pos[2]));
				energy+=fabs(Eng->GetVolt(2,pos[0],pos[1],pos[2]) * Eng->GetCurr(1,pos[0],pos[1],pos[2]));
			}
		}
	}
	return energy*0.5;
}

void  ProcessFields::SetSubSampling(unsigned int subSampleRate, int dir)
{
	if (dir>2) return;
	if (dir<0)
	{
		subSample[0]=subSampleRate;
		subSample[1]=subSampleRate;
		subSample[2]=subSampleRate;
	}
	else subSample[dir]=subSampleRate;
}

void ProcessFields::WriteVTKHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision)
{
	file << "# vtk DataFile Version 2.0" << endl;
	file << "Rectilinear Grid openEMS_ProcessFields" << endl;
	file << "ASCII" << endl;
	file << "DATASET RECTILINEAR_GRID " << endl;
	file << "DIMENSIONS " << numLines[0] << " " << numLines[1] << " " << numLines[2] << endl;
	file << "X_COORDINATES " << numLines[0] << " float" << endl;
	for (unsigned int i=0;i<numLines[0];++i)
		file << setprecision(precision) << discLines[0][i] << " ";
	file << endl;
	file << "Y_COORDINATES " << numLines[1] << " float" << endl;
	for (unsigned int i=0;i<numLines[1];++i)
		file << setprecision(precision) << discLines[1][i] << " ";
	file << endl;
	file << "Z_COORDINATES " << numLines[2] << " float" << endl;
	for (unsigned int i=0;i<numLines[2];++i)
		file << setprecision(precision) << discLines[2][i] << " ";
	file << endl << endl;
	file << "POINT_DATA " << numLines[0]*numLines[1]*numLines[2] << endl;
}

void ProcessFields::WriteVTKVectorArray(ofstream &file, string name, FDTD_FLOAT const* const* const* const* array, unsigned int const* numLines, unsigned int precision)
{
	file << "VECTORS " << name << " float " << endl;

	unsigned int pos[3];
	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				//in x
				file << setprecision(precision) << array[0][pos[0]][pos[1]][pos[2]] << " ";
				//in y
				file << setprecision(precision) << array[1][pos[0]][pos[1]][pos[2]] << " ";
				//in z
				file << setprecision(precision) << array[2][pos[0]][pos[1]][pos[2]] << endl;
			}
		}
	}
}


bool ProcessFields::DumpVectorArray2VTK(ofstream &file, string name, FDTD_FLOAT const* const* const* const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision)
{
	WriteVTKHeader(file, discLines, numLines, precision);
	WriteVTKVectorArray(file, name, array, numLines, precision);
	return true;
}

bool ProcessFields::DumpMultiVectorArray2VTK(ofstream &file, string names[], FDTD_FLOAT const* const* const* const* const* array, unsigned int numFields, double const* const* discLines, unsigned int const* numLines, unsigned int precision)
{
	WriteVTKHeader(file, discLines, numLines, precision);
	for (unsigned int n=0;n<numFields;++n)
	{
		WriteVTKVectorArray(file, names[n], array[n], numLines, precision);
		file << endl;
	}
	return true;
}

void ProcessFields::WriteVTKScalarArray(ofstream &file, string name, FDTD_FLOAT const* const* const* array, unsigned int const* numLines, unsigned int precision)
{
	file << "SCALARS " << name << " float " << 1 << endl;
	file << "LOOKUP_TABLE default" << endl;
	unsigned int pos[3];
	int count=0;
	for (pos[2]=0;pos[2]<numLines[2];++pos[2])
	{
		for (pos[1]=0;pos[1]<numLines[1];++pos[1])
		{
			for (pos[0]=0;pos[0]<numLines[0];++pos[0])
			{
				file << setprecision(precision) << array[pos[0]][pos[1]][pos[2]] << " ";
				++count;
				if (count%10==0)
					file << endl;
			}
		}
	}
}

bool ProcessFields::DumpScalarArray2VTK(ofstream &file, string name, FDTD_FLOAT const* const* const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision)
{
	WriteVTKHeader(file, discLines, numLines, precision);
	WriteVTKScalarArray(file, name, array, numLines, precision);
	return true;
}

bool ProcessFields::DumpMultiScalarArray2VTK(ofstream &file, string names[], FDTD_FLOAT const* const* const* const* array, unsigned int numFields, double const* const* discLines, unsigned int const* numLines, unsigned int precision)
{
	WriteVTKHeader(file, discLines, numLines);
	for (unsigned int n=0;n<numFields;++n)
	{
		WriteVTKScalarArray(file, names[n], array[n], numLines, precision);
		file << endl;
	}
	return true;
}

bool ProcessFields::DumpVectorArray2HDF5(string filename, string name, FDTD_FLOAT const* const* const* const* array, unsigned int const* numLines)
{
	const H5std_string FILE_NAME(filename);
	const H5std_string DATASET_NAME( name );

	H5::H5File file( FILE_NAME, H5F_ACC_RDWR );

	H5::Group group( file.openGroup( "/FieldData" ));

	hsize_t dimsf[4];              // dataset dimensions

	dimsf[0] = 3;
	dimsf[1] = numLines[2];
	dimsf[2] = numLines[1];
	dimsf[3] = numLines[0];

	H5::DataSpace dataspace( 4, dimsf );

	H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );
//	datatype.setOrder( H5T_ORDER_LE );
	H5::DataSet dataset = group.createDataSet( DATASET_NAME, datatype, dataspace );

	// I have not the slightest idea why this array-copy action is necessary...  but it's the only way hdf5 does what it is supposed to do anyway!!
	// at least it is save in case FDTD_FLOAT was defined as double...
	// why does hdf5 write the dimensions backwards??? or matlab???
	float hdf5array[3][numLines[2]][numLines[1]][numLines[0]];
	for (int n=0;n<3;++n)
	{
		for (unsigned int i=0;i<numLines[0];++i)
		{
			for (unsigned int j=0;j<numLines[1];++j)
			{
				for (unsigned int k=0;k<numLines[2];++k)
				{
					hdf5array[n][k][j][i] = array[n][i][j][k];
				}
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );
	return true;
}
