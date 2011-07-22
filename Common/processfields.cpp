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
#include "tools/vtk_file_io.h"
#include "processfields.h"
#include "FDTD/engine_interface_fdtd.h"

ProcessFields::ProcessFields(Engine_Interface_Base* eng_if) : Processing(eng_if)
{
	m_DumpType = E_FIELD_DUMP;
	// vtk-file is default
	m_fileType = VTK_FILETYPE;
	m_SampleType = NONE;
	m_Dump_File = NULL;
	SetPrecision(6);
	m_dualTime = false;

	// dump box should be always inside the snapped lines
	m_SnapMethod = 1;

	for (int n=0; n<3; ++n)
	{
		numLines[n]=0;
		posLines[n]=NULL;
		discLines[n]=NULL;
		subSample[n]=1;
		optResolution[n]=0;
	}
}

ProcessFields::~ProcessFields()
{
	delete m_Dump_File;
	m_Dump_File = NULL;
	for (int n=0; n<3; ++n)
	{
		delete[] posLines[n];
		posLines[n]=NULL;
		delete[] discLines[n];
		discLines[n]=NULL;
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
	case J_FIELD_DUMP:
		return "J-Field";
	case ROTH_FIELD_DUMP:
		return "RotH-Field";
	case SAR_LOCAL_DUMP:
		return "SAR-local";
	}
	return "unknown field";
}

void ProcessFields::InitProcess()
{
	if (Enabled==false) return;

	CalcMeshPos();

	if (m_fileType==VTK_FILETYPE)
	{
		delete m_Dump_File;
		m_Dump_File = new VTK_File_IO(m_filename,(int)m_Mesh_Type);

		#ifdef OUTPUT_IN_DRAWINGUNITS
		double discScaling = 1;
		#else
		double discScaling = Op->GetGridDelta();
		#endif
		m_Dump_File->SetMeshLines(discLines,numLines,discScaling);
		m_Dump_File->SetNativeDump(g_settings.NativeFieldDumps());
	}
}

void ProcessFields::SetDumpMode(Engine_Interface_Base::InterpolationType mode)
{
	m_Eng_Interface->SetInterpolationType(mode);
	if (mode==Engine_Interface_Base::CELL_INTERPOLATE)
		m_dualMesh=true;
	else if (mode==Engine_Interface_Base::NODE_INTERPOLATE)
		m_dualMesh=false;
	//else keep the preset/user defined case
}

void ProcessFields::DefineStartStopCoord(double* dstart, double* dstop)
{
	Processing::DefineStartStopCoord(dstart,dstop);

	// normalize order of start and stop
	for (int n=0; n<3; ++n)
	{
		if (start[n]>stop[n])
		{
			unsigned int help = start[n];
			start[n]=stop[n];
			stop[n]=help;
		}
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
		for (pos[0]=0; pos[0]<Op->GetNumberOfLines(0)-1; ++pos[0])
		{
			for (pos[1]=0; pos[1]<Op->GetNumberOfLines(1)-1; ++pos[1])
			{
				for (pos[2]=0; pos[2]<Op->GetNumberOfLines(2)-1; ++pos[2])
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
	m_SampleType = SUBSAMPLE;
}

void ProcessFields::SetOptResolution(double optRes, int dir)
{
	if (dir>2) return;
	if (dir<0)
	{
		optResolution[0]=optRes;
		optResolution[1]=optRes;
		optResolution[2]=optRes;
	}
	else optResolution[dir]=optRes;
	m_SampleType = OPT_RESOLUTION;
}

void ProcessFields::CalcMeshPos()
{
	if ((m_SampleType==SUBSAMPLE) || (m_SampleType==NONE))
	{
		vector<unsigned int> tmp_pos;

		for (int n=0; n<3; ++n)
		{
			// construct new discLines
			tmp_pos.clear();
			for (unsigned int i=start[n]; i<=stop[n]; i+=subSample[n])
				tmp_pos.push_back(i);

			numLines[n] = tmp_pos.size();
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			delete[] posLines[n];
			posLines[n] = new unsigned int[numLines[n]];
			for (unsigned int i=0; i<numLines[n]; ++i)
			{
				posLines[n][i] = tmp_pos.at(i);
				discLines[n][i] = Op->GetDiscLine(n,tmp_pos.at(i),m_dualMesh);
			}
		}
	}
	if ((m_SampleType==OPT_RESOLUTION))
	{
		vector<unsigned int> tmp_pos;
		double oldPos=0;
		for (int n=0; n<3; ++n)
		{
			// construct new discLines
			tmp_pos.clear();
			tmp_pos.push_back(start[n]);
			oldPos=Op->GetDiscLine(n,start[n],m_dualMesh);
			if (stop[n]==0)
				tmp_pos.push_back(stop[n]);
			else
				for (unsigned int i=start[n]+1; i<=stop[n]-1; ++i)
				{
					if ( (Op->GetDiscLine(n,i+1,m_dualMesh)-oldPos) >= optResolution[n])
					{
						tmp_pos.push_back(i);
						oldPos=Op->GetDiscLine(n,i,m_dualMesh);
					}
				}
			if (start[n]!=stop[n])
				tmp_pos.push_back(stop[n]);
			numLines[n] = tmp_pos.size();
			delete[] discLines[n];
			discLines[n] = new double[numLines[n]];
			delete[] posLines[n];
			posLines[n] = new unsigned int[numLines[n]];
			for (unsigned int i=0; i<numLines[n]; ++i)
			{
				posLines[n][i] = tmp_pos.at(i);
				discLines[n][i] = Op->GetDiscLine(n,tmp_pos.at(i),m_dualMesh);
			}
		}
	}
}

bool ProcessFields::WriteMesh2HDF5(string filename, string groupName, unsigned int const* numLines, double const* const* discLines, MeshType meshT, double discLines_scaling)
{
	H5::H5File file( filename, H5F_ACC_RDWR );

	H5::Group hdf_group( file.openGroup( groupName ));

	string names[] = {"x","y","z"};
	if (meshT==CYLINDRICAL_MESH)
	{
		names[0]="rho";
		names[1]="alpha";
	}

	H5::Group* group = new H5::Group( hdf_group.createGroup( "/Mesh" ));
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
			if ((meshT==CYLINDRICAL_MESH) && (n==1)) //check for alpha-direction
				array[i] = discLines[n][i];
			else
				array[i] = discLines[n][i] * discLines_scaling;
#endif
		}
		//write to dataset
		dataset.write( array, H5::PredType::NATIVE_FLOAT );
		delete[] array;
	}
	delete group;

	return true;
}

bool ProcessFields::DumpVectorArray2HDF5(string filename, string groupName, string name, FDTD_FLOAT const* const* const* const* array, unsigned int const* numLines, float time)
{
	const H5std_string FILE_NAME(filename);
	const H5std_string DATASET_NAME( name );

	H5::H5File file( FILE_NAME, H5F_ACC_RDWR );

	H5::Group group( file.openGroup( groupName ));

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
	unsigned long pos = 0;
	float *hdf5array = new float[3*numLines[0]*numLines[1]*numLines[2]];
	for (int n=0; n<3; ++n)
	{
		for (unsigned int k=0; k<numLines[2]; ++k)
		{
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				for (unsigned int i=0; i<numLines[0]; ++i)
				{
					hdf5array[pos++] = array[n][i][j][k];
				}
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );
	delete[] hdf5array;
	return true;
}

bool ProcessFields:: DumpScalarArray2HDF5(string filename, string groupName, string name, FDTD_FLOAT const* const* const* array, unsigned int const* numLines, string attr_name, float attr_value)
{
	const H5std_string FILE_NAME(filename);
	const H5std_string DATASET_NAME( name );

	H5::H5File file( FILE_NAME, H5F_ACC_RDWR );

	H5::Group group( file.openGroup( groupName ));

	hsize_t dimsf[3];              // dataset dimensions

	dimsf[0] = numLines[2];
	dimsf[1] = numLines[1];
	dimsf[2] = numLines[0];

	H5::DataSpace dataspace( 3, dimsf );

	H5::FloatType datatype( H5::PredType::NATIVE_FLOAT );
//	datatype.setOrder( H5T_ORDER_LE );
	H5::DataSet dataset = group.createDataSet( DATASET_NAME, datatype, dataspace );

	if (!attr_name.empty())
	{
		hsize_t t_dimsf[] = {1};
		H5::DataSpace t_dataspace( 1, t_dimsf );
		H5::Attribute attr = dataset.createAttribute(attr_name,H5::PredType::NATIVE_FLOAT,t_dataspace);
		attr.write( H5::PredType::NATIVE_FLOAT , &attr_value);
	}

	// I have not the slightest idea why this array-copy action is necessary...  but it's the only way hdf5 does what it is supposed to do anyway!!
	// at least it is save in case FDTD_FLOAT was defined as double...
	// why does hdf5 write the dimensions backwards??? or matlab???
	unsigned long pos = 0;
	float *hdf5array = new float[numLines[0]*numLines[1]*numLines[2]];
	for (unsigned int k=0; k<numLines[2]; ++k)
	{
		for (unsigned int j=0; j<numLines[1]; ++j)
		{
			for (unsigned int i=0; i<numLines[0]; ++i)
			{
				hdf5array[pos++] = array[i][j][k];
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );
	delete[] hdf5array;
	return true;
}

bool ProcessFields::DumpVectorArray2HDF5(string filename, string groupName, string name, std::complex<float> const* const* const* const* array, unsigned int const* numLines, float weight, float frequency)
{
	const H5std_string FILE_NAME(filename);
	const H5std_string DATASET_NAME_RE( name + "_real");
	const H5std_string DATASET_NAME_IM( name + "_imag");

	H5::H5File file( FILE_NAME, H5F_ACC_RDWR );

	H5::Group group( file.openGroup( groupName ));

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
	unsigned long pos = 0;
	float *hdf5array = new float[3*numLines[0]*numLines[1]*numLines[2]];
	for (int n=0; n<3; ++n)
	{
		for (unsigned int k=0; k<numLines[2]; ++k)
		{
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				for (unsigned int i=0; i<numLines[0]; ++i)
				{
					hdf5array[pos++] = array[n][i][j][k].real() * weight;
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
	pos=0;
	for (int n=0; n<3; ++n)
	{
		for (unsigned int k=0; k<numLines[2]; ++k)
		{
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				for (unsigned int i=0; i<numLines[0]; ++i)
				{
					hdf5array[pos++] = array[n][i][j][k].imag() * weight;
				}
			}
		}
	}
	dataset.write( hdf5array, H5::PredType::NATIVE_FLOAT );
	delete[] hdf5array;
	return true;
}

FDTD_FLOAT**** ProcessFields::CalcField()
{
	unsigned int pos[3];
	double out[3];
	//create array
	FDTD_FLOAT**** field = Create_N_3DArray<FDTD_FLOAT>(numLines);
	switch (m_DumpType)
	{
	case E_FIELD_DUMP:
	case SAR_LOCAL_DUMP:
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			pos[0]=posLines[0][i];
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				pos[1]=posLines[1][j];
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					pos[2]=posLines[2][k];

					m_Eng_Interface->GetEField(pos,out);
					field[0][i][j][k] = out[0];
					field[1][i][j][k] = out[1];
					field[2][i][j][k] = out[2];
				}
			}
		}
		return field;
	case H_FIELD_DUMP:
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			pos[0]=posLines[0][i];
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				pos[1]=posLines[1][j];
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					pos[2]=posLines[2][k];

					m_Eng_Interface->GetHField(pos,out);
					field[0][i][j][k] = out[0];
					field[1][i][j][k] = out[1];
					field[2][i][j][k] = out[2];
				}
			}
		}
		return field;
	case J_FIELD_DUMP:
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			pos[0]=posLines[0][i];
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				pos[1]=posLines[1][j];
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					pos[2]=posLines[2][k];

					m_Eng_Interface->GetJField(pos,out);
					field[0][i][j][k] = out[0];
					field[1][i][j][k] = out[1];
					field[2][i][j][k] = out[2];
				}
			}
		}
		return field;
	case ROTH_FIELD_DUMP:
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			pos[0]=posLines[0][i];
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				pos[1]=posLines[1][j];
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					pos[2]=posLines[2][k];

					m_Eng_Interface->GetRotHField(pos,out);
					field[0][i][j][k] = out[0];
					field[1][i][j][k] = out[1];
					field[2][i][j][k] = out[2];
				}
			}
		}
		return field;
	}
	cerr << "ProcessFields::CalcField(): Error, unknown dump type..." << endl;
	return field;
}

