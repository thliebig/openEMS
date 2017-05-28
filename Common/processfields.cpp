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
#include "tools/global.h"
#include "tools/vtk_file_writer.h"
#include "tools/hdf5_file_writer.h"
#include "processfields.h"
#include "FDTD/engine_interface_fdtd.h"

ProcessFields::ProcessFields(Engine_Interface_Base* eng_if) : Processing(eng_if)
{
	m_DumpType = E_FIELD_DUMP;
	// vtk-file is default
	m_fileType = VTK_FILETYPE;
	m_SampleType = NONE;
	m_Vtk_Dump_File = NULL;
	m_HDF5_Dump_File = NULL;
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
	delete m_Vtk_Dump_File;
	m_Vtk_Dump_File = NULL;
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
	case D_FIELD_DUMP:
		return "D-Field";
	case B_FIELD_DUMP:
		return "B-Field";
	case SAR_LOCAL_DUMP:
		return "SAR-local";
	case SAR_1G_DUMP:
		return "SAR_1g";
	case SAR_10G_DUMP:
		return "SAR_10g";
	case SAR_RAW_DATA:
		return "SAR_raw_data";
	}
	return "unknown field";
}

bool ProcessFields::NeedConductivity() const
{
	switch (m_DumpType)
	{
	case J_FIELD_DUMP:
		return true;
	default:
		return false;
	}
	return false;
}

bool ProcessFields::NeedPermittivity() const
{
	switch (m_DumpType)
	{
	case D_FIELD_DUMP:
		return true;
	default:
		return false;
	}
	return false;
}

bool ProcessFields::NeedPermeability() const
{
	switch (m_DumpType)
	{
	case B_FIELD_DUMP:
		return true;
	default:
		return false;
	}
	return false;
}

void ProcessFields::InitProcess()
{
	if (Enabled==false) return;

	CalcMeshPos();

	if (m_fileType==VTK_FILETYPE)
	{
		delete m_Vtk_Dump_File;
		m_Vtk_Dump_File = new VTK_File_Writer(m_filename,(int)m_Mesh_Type);

		#ifdef OUTPUT_IN_DRAWINGUNITS
		double discScaling = 1;
		#else
		double discScaling = Op->GetGridDelta();
		#endif
		m_Vtk_Dump_File->SetMeshLines(discLines,numLines,discScaling);
		m_Vtk_Dump_File->SetNativeDump(g_settings.NativeFieldDumps());
	}
	if (m_fileType==HDF5_FILETYPE)
	{
		delete m_HDF5_Dump_File;
		m_HDF5_Dump_File = new HDF5_File_Writer(m_filename+".h5");

		#ifdef OUTPUT_IN_DRAWINGUNITS
		double discScaling = 1;
		#else
		double discScaling = Op->GetGridDelta();
		#endif
		m_HDF5_Dump_File->WriteRectMesh(numLines,discLines,(int)m_Mesh_Type,discScaling);

		m_HDF5_Dump_File->WriteAtrribute("/","openEMS_HDF5_version",0.2);
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

double ProcessFields::CalcTotalEnergyEstimate() const
{
	return m_Eng_Interface->CalcFastEnergy();
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

FDTD_FLOAT**** ProcessFields::CalcField()
{
	unsigned int pos[3];
	double out[3];
	//create array
	FDTD_FLOAT**** field = Create_N_3DArray<FDTD_FLOAT>(numLines);
	switch (m_DumpType)
	{
	case E_FIELD_DUMP:
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
	case D_FIELD_DUMP:
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			pos[0]=posLines[0][i];
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				pos[1]=posLines[1][j];
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					pos[2]=posLines[2][k];

					m_Eng_Interface->GetDField(pos,out);
					field[0][i][j][k] = out[0];
					field[1][i][j][k] = out[1];
					field[2][i][j][k] = out[2];
				}
			}
		}
		return field;
	case B_FIELD_DUMP:
		for (unsigned int i=0; i<numLines[0]; ++i)
		{
			pos[0]=posLines[0][i];
			for (unsigned int j=0; j<numLines[1]; ++j)
			{
				pos[1]=posLines[1][j];
				for (unsigned int k=0; k<numLines[2]; ++k)
				{
					pos[2]=posLines[2][k];

					m_Eng_Interface->GetBField(pos,out);
					field[0][i][j][k] = out[0];
					field[1][i][j][k] = out[1];
					field[2][i][j][k] = out[2];
				}
			}
		}
		return field;
	default:
		cerr << "ProcessFields::CalcField(): Error, unknown dump type..." << endl;
		return field;
	}
}

