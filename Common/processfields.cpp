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

#include <iomanip>
#include <H5Cpp.h>
#include "tools/global.h"
#include "processfields.h"
#include "FDTD/engine_interface_fdtd.h"

ProcessFields::ProcessFields(Engine_Interface_Base* eng_if) : Processing(eng_if)
{
	m_DumpType = E_FIELD_DUMP;
	// vtk-file is default
	m_fileType = VTK_FILETYPE;
	SetSubSampling(1);
	SetPrecision(6);
	m_dualTime = false;

	for (int n=0; n<3; ++n)
	{
		numLines[n]=0;
		discLines[n]=NULL;
	}
}

ProcessFields::~ProcessFields()
{
	for (int n=0; n<3; ++n)
	{
		delete[] discLines[n];
		discLines[n]=NULL;
	}
}

void ProcessFields::InitProcess()
{
	if (Enabled==false) return;
	//get the correct direction names for all coordinate systems
	string names[] = {Op->GetDirName(0),Op->GetDirName(1),Op->GetDirName(2)};
	if (m_fileType==HDF5_FILETYPE)
	{
		m_filename+= ".h5";

		H5::H5File* file = new H5::H5File( m_filename , H5F_ACC_TRUNC );

		H5::Group* group = new H5::Group( file->createGroup( "/Mesh" ));
		for (int n=0; n<3; ++n)
		{
			hsize_t dimsf[1];              // dataset dimensions
			dimsf[0] = numLines[n];
			H5::DataSpace dataspace( 1, dimsf );
			H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );
			H5::DataSet dataset = group->createDataSet( names[n].c_str(), datatype, dataspace );
			//convert to float...
			float* array = new float[numLines[n]];
			for (unsigned int i=0; i<numLines[n]; ++i)
			{
#ifdef OUTPUT_IN_DRAWINGUNITS
				array[i] = Lines[n][i];
#else
				if ((m_Mesh_Type==CYLINDRICAL_MESH) && (n==1)) //check for alpha-direction
					array[i] = discLines[n][i];
				else
					array[i] = discLines[n][i] * Op->GetGridDelta();
#endif
			}
			//write to dataset
			dataset.write( array, H5::PredType::NATIVE_FLOAT );
		}
		delete group;

		group = new H5::Group( file->createGroup( "/FieldData" ));
		delete group;
		group = new H5::Group( file->createGroup( "/FieldData/FD" ));
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

	// determine mesh type
	bool dualMesh = false;
	if (m_DumpType == H_FIELD_DUMP)
		dualMesh = true;

	Engine_Interface_Base::InterpolationType m_DumpMode = m_Eng_Interface->GetInterpolationType();
	if (m_DumpMode==Engine_Interface_Base::NO_INTERPOLATION)
	{
		if (!Op->SnapToMesh(dstart,start,dualMesh))
			cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (!Op->SnapToMesh(dstop,stop,dualMesh))
			cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;

		for (int n=0; n<3; ++n)
		{
			// normalize order of start and stop
			if (start[n]>stop[n])
			{
				unsigned int help = start[n];
				start[n]=stop[n];
				stop[n]=help;
			}

			// construct new discLines
			lines.clear();
			for (unsigned int i=start[n]; i<=stop[n]; i+=subSample[n])
			{
				lines.push_back(Op->GetDiscLine(n,i,dualMesh));
			}
			numLines[n] = lines.size();
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			for (unsigned int i=0; i<numLines[n]; ++i)
				discLines[n][i] = lines.at(i);
		}
	}
	else if (m_DumpMode==Engine_Interface_Base::NODE_INTERPOLATE)
	{
		if (Op->SnapToMesh(dstart,start)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (Op->SnapToMesh(dstop,stop)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;

		//create mesh
		for (int n=0; n<3; ++n)
		{
			if (start[n]>stop[n])
			{
				unsigned int help = start[n];
				start[n]=stop[n];
				stop[n]=help;
			}
//			if (stop[n]==Op->GetNumberOfLines(n)-1)
//				--stop[n];
//			cerr << "start " << start[n] << "stop " << stop[n];
			lines.clear();
			for (unsigned int i=start[n]; i<=stop[n]; i+=subSample[n])
			{
				lines.push_back(Op->GetDiscLine(n,i));//0.5*(Op->discLines[n][i+1] +  Op->discLines[n][i]));
			}
			numLines[n] = lines.size();
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			for (unsigned int i=0; i<numLines[n]; ++i)
				discLines[n][i] = lines.at(i);
		}
	}
	else if (m_DumpMode==Engine_Interface_Base::CELL_INTERPOLATE)
	{
		if (Op->SnapToMesh(dstart,start,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check start value!!" << endl;
		if (Op->SnapToMesh(dstop,stop,true)==false) cerr << "ProcessFields::DefineStartStopCoord: Warning: Snapping problem, check stop value!!" << endl;

		//create dual mesh
		for (int n=0; n<3; ++n)
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
			for (unsigned int i=start[n]; i<stop[n]; i+=subSample[n])
			{
				lines.push_back(Op->GetDiscLine(n,i,true));//0.5*(Op->discLines[n][i+1] +  Op->discLines[n][i]));
			}
			numLines[n] = lines.size();
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			for (unsigned int i=0; i<numLines[n]; ++i)
				discLines[n][i] = lines.at(i);
		}
	}

	if (g_settings.showProbeDiscretization())
	{
		// FIXME the information E-Field / H-Field and therefore which mesh to use is missing
		bool dualMesh = false;
		cerr << m_filename << ": snapped coords: (" << Op->GetDiscLine( 0, start[0], dualMesh ) << ","
		     << Op->GetDiscLine( 1, start[1], dualMesh ) << "," << Op->GetDiscLine( 2, start[2], dualMesh ) << ") -> ("
		     << Op->GetDiscLine( 0, stop[0], dualMesh ) << ","<< Op->GetDiscLine( 1, stop[1], dualMesh ) << ","
		     << Op->GetDiscLine( 2, stop[2], dualMesh ) << ")";
		cerr << "   [" << start[0] << "," << start[1] << "," << start[2] << "] -> ["
		     << stop[0] << "," << stop[1] << "," << stop[2] << "]" << endl;
	}

}

double ProcessFields::CalcTotalEnergy() const
{
	double energy=0.0;

	Engine_Interface_FDTD* EI_FDTD = dynamic_cast<Engine_Interface_FDTD*>(m_Eng_Interface);

	if (EI_FDTD)
	{
		const Engine* Eng = EI_FDTD->GetFDTDEngine();

		unsigned int pos[3];
		for (pos[0]=0; pos[0]<Op->GetNumberOfLines(0); ++pos[0])
		{
			for (pos[1]=0; pos[1]<Op->GetNumberOfLines(1); ++pos[1])
			{
				for (pos[2]=0; pos[2]<Op->GetNumberOfLines(2); ++pos[2])
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
	}
	return energy*0.5;
}

void ProcessFields::SetSubSampling(unsigned int subSampleRate, int dir)
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

void ProcessFields::WriteVTKHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, MeshType meshT, double discLines_scaling)
{
	if (meshT==CARTESIAN_MESH)
		WriteVTKCartesianGridHeader(file, discLines, numLines, precision, header_info, discLines_scaling);
	else if (meshT==CYLINDRICAL_MESH)
		WriteVTKCylindricalGridHeader(file, discLines, numLines, precision, header_info, discLines_scaling);
	else
		cerr << "ProcessFields::WriteVTKHeader: Warning: unknown mesh type, skipping header -> file will be invalid..." << endl;
}

void ProcessFields::WriteVTKCartesianGridHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, double discLines_scaling)
{
	file << "# vtk DataFile Version 2.0" << endl;
	file << "Rectilinear Grid openEMS_ProcessFields";
	if (!header_info.empty())
		file << " " << header_info;
	file << endl;
	file << "ASCII" << endl;
	file << "DATASET RECTILINEAR_GRID " << endl;
	file << "DIMENSIONS " << numLines[0] << " " << numLines[1] << " " << numLines[2] << endl;
	file << "X_COORDINATES " << numLines[0] << " " << __VTK_DATA_TYPE__ << endl;
	for (unsigned int i=0; i<numLines[0]; ++i)
		file << setprecision(precision) << discLines[0][i] * discLines_scaling << " ";
	file << endl;
	file << "Y_COORDINATES " << numLines[1] << " " << __VTK_DATA_TYPE__ << endl;
	for (unsigned int i=0; i<numLines[1]; ++i)
		file << setprecision(precision) << discLines[1][i] * discLines_scaling << " ";
	file << endl;
	file << "Z_COORDINATES " << numLines[2] << " " << __VTK_DATA_TYPE__ << endl;
	for (unsigned int i=0; i<numLines[2]; ++i)
		file << setprecision(precision) << discLines[2][i] * discLines_scaling << " ";
	file << endl << endl;
	file << "POINT_DATA " << numLines[0]*numLines[1]*numLines[2] << endl;
}

void ProcessFields::WriteVTKCylindricalGridHeader(ofstream &file, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, double discLines_scaling)
{
	file << "# vtk DataFile Version 3.0" << endl;
	file << "Structured Grid from openEMS_ProcessFields";
	if (!header_info.empty())
		file << " " << header_info;
	file << endl;
	file << "ASCII" << endl;
	file << "DATASET STRUCTURED_GRID " << endl;
	file << "DIMENSIONS " << numLines[0] << " " << numLines[1] << " " << numLines[2] << endl;
	file << "POINTS " << numLines[0]*numLines[1]*numLines[2] << " " << __VTK_DATA_TYPE__ << endl;
	for (unsigned int k=0; k<numLines[2]; ++k)
		for (unsigned int j=0; j<numLines[1]; ++j)
			for (unsigned int i=0; i<numLines[0]; ++i)
			{
				file << setprecision(precision) << discLines[0][i] * cos(discLines[1][j]) * discLines_scaling << " "
				<< discLines[0][i] * sin(discLines[1][j]) * discLines_scaling << " "
				<< discLines[2][k] * discLines_scaling << endl;
			}
	file << endl;
	file << endl << endl;
	file << "POINT_DATA " << numLines[0]*numLines[1]*numLines[2] << endl;
}


void ProcessFields::WriteVTKVectorArray(ofstream &file, string name, FDTD_FLOAT const* const* const* const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision, MeshType meshT)
{
	file << "VECTORS " << name << " " << __VTK_DATA_TYPE__ << endl;

	if (g_settings.NativeFieldDumps())
		meshT = CARTESIAN_MESH; //dump field components as they are...

	unsigned int pos[3];
	for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			double cos_a = cos(discLines[1][pos[1]]); //needed only for CYLINDRICAL_MESH
			double sin_a = sin(discLines[1][pos[1]]); //needed only for CYLINDRICAL_MESH
			for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			{
				switch (meshT)
				{
				case CARTESIAN_MESH:
					UNUSED(discLines); //disclines not needed for the original cartesian mesh
					//in x
					file << setprecision(precision) << array[0][pos[0]][pos[1]][pos[2]] << " ";
					//in y
					file << setprecision(precision) << array[1][pos[0]][pos[1]][pos[2]] << " ";
					//in z
					file << setprecision(precision) << array[2][pos[0]][pos[1]][pos[2]] << endl;
					break;
				case CYLINDRICAL_MESH:
					//in x : F_x = F_r * cos(a) - F_a * sin(a);
					file << setprecision(precision) << array[0][pos[0]][pos[1]][pos[2]] * cos_a - array[1][pos[0]][pos[1]][pos[2]] * sin_a << " ";
					//in y : F_y = F_r * sin(a) + F_a * cos(a);
					file << setprecision(precision) << array[0][pos[0]][pos[1]][pos[2]] * sin_a + array[1][pos[0]][pos[1]][pos[2]] * cos_a << " ";
					//in z
					file << setprecision(precision) << array[2][pos[0]][pos[1]][pos[2]] << endl;
					break;
				}
			}
		}
	}
}


bool ProcessFields::DumpVectorArray2VTK(ofstream &file, string name, FDTD_FLOAT const* const* const* const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, MeshType meshT, double discLines_scaling)
{
	WriteVTKHeader(file, discLines, numLines, precision, header_info, meshT, discLines_scaling);
	WriteVTKVectorArray(file, name, array, discLines, numLines, precision, meshT);
	return true;
}

bool ProcessFields::DumpMultiVectorArray2VTK(ofstream &file, string names[], FDTD_FLOAT const* const* const* const* const* array, unsigned int numFields, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, MeshType meshT, double discLines_scaling)
{
	WriteVTKHeader(file, discLines, numLines, precision, header_info, meshT, discLines_scaling);
	for (unsigned int n=0; n<numFields; ++n)
	{
		WriteVTKVectorArray(file, names[n], array[n], discLines, numLines, precision, meshT);
		file << endl;
	}
	return true;
}

void ProcessFields::WriteVTKScalarArray(ofstream &file, string name, FDTD_FLOAT const* const* const* array, unsigned int const* numLines, unsigned int precision)
{
	file << "SCALARS " << name << " " << __VTK_DATA_TYPE__ << 1 << endl;
	file << "LOOKUP_TABLE default" << endl;
	unsigned int pos[3];
	int count=0;
	for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
	{
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			{
				file << setprecision(precision) << array[pos[0]][pos[1]][pos[2]] << " ";
				++count;
				if (count%10==0)
					file << endl;
			}
		}
	}
}

bool ProcessFields::DumpScalarArray2VTK(ofstream &file, string name, FDTD_FLOAT const* const* const* array, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, MeshType meshT, double discLines_scaling)
{
	WriteVTKHeader(file, discLines, numLines, precision, header_info, meshT, discLines_scaling);
	WriteVTKScalarArray(file, name, array, numLines, precision);
	return true;
}

bool ProcessFields::DumpMultiScalarArray2VTK(ofstream &file, string names[], FDTD_FLOAT const* const* const* const* array, unsigned int numFields, double const* const* discLines, unsigned int const* numLines, unsigned int precision, string header_info, MeshType meshT, double discLines_scaling)
{
	WriteVTKHeader(file, discLines, numLines, precision, header_info, meshT, discLines_scaling);
	for (unsigned int n=0; n<numFields; ++n)
	{
		WriteVTKScalarArray(file, names[n], array[n], numLines, precision);
		file << endl;
	}
	return true;
}

bool ProcessFields::DumpVectorArray2HDF5(string filename, string name, FDTD_FLOAT const* const* const* const* array, unsigned int const* numLines, float time)
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

	hsize_t t_dimsf[] = {1};
	H5::DataSpace t_dataspace( 1, t_dimsf );
	H5::Attribute attr = dataset.createAttribute("time",H5::PredType::NATIVE_FLOAT,t_dataspace);
	attr.write( H5::PredType::NATIVE_FLOAT , &time);

	// I have not the slightest idea why this array-copy action is necessary...  but it's the only way hdf5 does what it is supposed to do anyway!!
	// at least it is save in case FDTD_FLOAT was defined as double...
	// why does hdf5 write the dimensions backwards??? or matlab???
	float hdf5array[3][numLines[2]][numLines[1]][numLines[0]];
	for (int n=0; n<3; ++n)
	{
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					hdf5array[n][k][j][i] = array[n][i][j][k];
				}
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );
	return true;
}

bool ProcessFields::DumpVectorArray2HDF5(string filename, string name, std::complex<float> const* const* const* const* array, unsigned int const* numLines, float weight, float frequency)
{
	const H5std_string FILE_NAME(filename);
	const H5std_string DATASET_NAME_RE( name + "_real");
	const H5std_string DATASET_NAME_IM( name + "_imag");

	H5::H5File file( FILE_NAME, H5F_ACC_RDWR );

	H5::Group group( file.openGroup( "/FieldData/FD" ));

	hsize_t t_dimsf[] = {1};
	H5::DataSpace t_dataspace( 1, t_dimsf );

	hsize_t dimsf[4];              // dataset dimensions
	dimsf[0] = 3;
	dimsf[1] = numLines[2];
	dimsf[2] = numLines[1];
	dimsf[3] = numLines[0];

	H5::DataSpace dataspace( 4, dimsf );
	H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );

	//create and write real part
	H5::DataSet dataset = group.createDataSet( DATASET_NAME_RE, datatype, dataspace );
	H5::Attribute attr = dataset.createAttribute("frequency",H5::PredType::NATIVE_FLOAT,t_dataspace);
	attr.write( H5::PredType::NATIVE_FLOAT , &frequency);
	// I have not the slightest idea why this array-copy action is necessary...  but it's the only way hdf5 does what it is supposed to do anyway!!
	// at least it is save in case FDTD_FLOAT was defined as double...
	// why does hdf5 write the dimensions backwards??? or matlab???
	float hdf5array[3][numLines[2]][numLines[1]][numLines[0]];
	for (int n=0; n<3; ++n)
	{
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					hdf5array[n][k][j][i] = array[n][i][j][k].real() * weight;
				}
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );

	//create and write imaginary part
	dataset = group.createDataSet( DATASET_NAME_IM, datatype, dataspace );
	attr = dataset.createAttribute("frequency",H5::PredType::NATIVE_FLOAT,t_dataspace);
	attr.write( H5::PredType::NATIVE_FLOAT , &frequency);
	// I have not the slightest idea why this array-copy action is necessary...  but it's the only way hdf5 does what it is supposed to do anyway!!
	// at least it is save in case FDTD_FLOAT was defined as double...
	// why does hdf5 write the dimensions backwards??? or matlab???
	for (int n=0; n<3; ++n)
	{
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					hdf5array[n][k][j][i] = array[n][i][j][k].imag() * weight;
				}
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );

	return true;
}

FDTD_FLOAT**** ProcessFields::CalcField()
{
	unsigned int pos[3];
	unsigned int OpPos[3];
	double out[3];
	//create array
	FDTD_FLOAT**** field = Create_N_3DArray<FDTD_FLOAT>(numLines);
	if (m_DumpType==E_FIELD_DUMP)
	{
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					m_Eng_Interface->GetEField(OpPos,out);
					field[0][pos[0]][pos[1]][pos[2]] = out[0];
					field[1][pos[0]][pos[1]][pos[2]] = out[1];
					field[2][pos[0]][pos[1]][pos[2]] = out[2];
				}
			}
		}
		return field;
	}

	if (m_DumpType==H_FIELD_DUMP)
	{
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					m_Eng_Interface->GetHField(OpPos,out);
					field[0][pos[0]][pos[1]][pos[2]] = out[0];
					field[1][pos[0]][pos[1]][pos[2]] = out[1];
					field[2][pos[0]][pos[1]][pos[2]] = out[2];
				}
			}
		}
		return field;
	}

	cerr << "ProcessFields::CalcField(): Error, unknown dump type..." << endl;
	return field;
}

