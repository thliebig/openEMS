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

#include "processfields_td.h"
#include "Common/operator_base.h"
#include <iomanip>
#include <sstream>
#include <string>

ProcessFieldsTD::ProcessFieldsTD(Engine_Interface_Base* eng_if) : ProcessFields(eng_if)
{
	pad_length = 8;
}

ProcessFieldsTD::~ProcessFieldsTD()
{
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (filePattern.empty()) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	string filename;

	if (m_fileType==VTK_FILETYPE)
	{
		stringstream ss;
		ss << filePattern << std::setw( pad_length ) << std::setfill( '0' ) << m_Eng_Interface->GetNumberOfTimesteps() << ".vtk";
		filename = ss.str();
	}
	else
		filename = m_filename;

#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = Op->GetGridDelta();
#endif

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
	}

	if (m_fileType==VTK_FILETYPE)
	{
		ofstream file(filename.c_str());
		if (file.is_open()==false)
		{
			cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;
		};
		DumpVectorArray2VTK(file,GetFieldNameByType(m_DumpType),field,discLines,numLines,m_precision,string("Interpolation: ")+m_Eng_Interface->GetInterpolationTypeString(), m_Mesh_Type, discLines_scaling);
		file.close();
	}
	else if (m_fileType==HDF5_FILETYPE)
	{
		stringstream ss;
		ss << std::setw( pad_length ) << std::setfill( '0' ) << m_Eng_Interface->GetNumberOfTimesteps();
		DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),field,numLines,(0.5+m_Eng_Interface->GetNumberOfTimesteps())*Op->GetTimestep());
	}
	else
		cerr << "ProcessFieldsTD::Process: unknown File-Type" << endl;

	Delete_N_3DArray(field,numLines);


	return GetNextInterval();
}
