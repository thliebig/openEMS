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
#include "tools/vtk_file_io.h"
#include <H5Cpp.h>
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

void ProcessFieldsTD::InitProcess()
{
	if (Enabled==false) return;
	
	ProcessFields::InitProcess();

	if (m_Dump_File)
		m_Dump_File->SetHeader(string("openEMS TD Field Dump -- Interpolation: ")+m_Eng_Interface->GetInterpolationTypeString());

	#ifdef OUTPUT_IN_DRAWINGUNITS
	double discScaling = 1;
	#else
	double discScaling = Op->GetGridDelta();
	#endif

	if (m_fileType==HDF5_FILETYPE)
	{
		//create hdf5 file & necessary groups
		m_filename+= ".h5";
		H5::H5File* file = new H5::H5File( m_filename, H5F_ACC_TRUNC );
		H5::Group* group = new H5::Group( file->createGroup( "/FieldData" ));
		delete group;
		group = new H5::Group( file->createGroup( "/FieldData/TD" ));
		delete group;
		delete file;

		//write mesh information in main root-group
		ProcessFields::WriteMesh2HDF5(m_filename,"/",numLines,discLines,m_Mesh_Type, discScaling);
	}
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	string filename = m_filename;

	float**** field = CalcField();

	if (m_fileType==VTK_FILETYPE)
	{
		m_Dump_File->SetTimestep(m_Eng_Interface->GetNumberOfTimesteps());
		m_Dump_File->ClearAllFields();
		m_Dump_File->AddVectorField(GetFieldNameByType(m_DumpType),field,numLines);
		if (m_Dump_File->Write()==false)
			cerr << "ProcessFieldsTD::Process: can't dump to file... abort! " << endl;
	}
	else if (m_fileType==HDF5_FILETYPE)
	{
		stringstream ss;
		ss << std::setw( pad_length ) << std::setfill( '0' ) << m_Eng_Interface->GetNumberOfTimesteps();
		DumpVectorArray2HDF5(filename.c_str(), "/FieldData/TD", string( ss.str() ), field, numLines, m_Eng_Interface->GetTime(m_dualTime));
	}
	else
		cerr << "ProcessFieldsTD::Process: unknown File-Type" << endl;

	Delete_N_3DArray<FDTD_FLOAT>(field,numLines);

	return GetNextInterval();
}
