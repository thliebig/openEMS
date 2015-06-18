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
#include "tools/vtk_file_writer.h"
#include "tools/hdf5_file_writer.h"
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

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

	if (m_Vtk_Dump_File)
		m_Vtk_Dump_File->SetHeader(string("openEMS TD Field Dump -- Interpolation: ")+m_Eng_Interface->GetInterpolationTypeString());

	if (m_HDF5_Dump_File)
		m_HDF5_Dump_File->SetCurrentGroup("/FieldData/TD");
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	string filename = m_filename;

	float**** field = CalcField();
	bool success = true;

	if (m_fileType==VTK_FILETYPE)
	{
		m_Vtk_Dump_File->SetTimestep(m_Eng_Interface->GetNumberOfTimesteps());
		m_Vtk_Dump_File->ClearAllFields();
		m_Vtk_Dump_File->AddVectorField(GetFieldNameByType(m_DumpType),field);
		success &= m_Vtk_Dump_File->Write();
	}
	else if (m_fileType==HDF5_FILETYPE)
	{
		stringstream ss;
		ss << std::setw( pad_length ) << std::setfill( '0' ) << m_Eng_Interface->GetNumberOfTimesteps();
		size_t datasize[]={numLines[0],numLines[1],numLines[2]};
		success &= m_HDF5_Dump_File->WriteVectorField(ss.str(), field, datasize);
		float time[1] = {(float)m_Eng_Interface->GetTime(m_dualTime)};
		success &= m_HDF5_Dump_File->WriteAtrribute("/FieldData/TD/"+ss.str(),"time",time,1);
	}
	else
	{
		success = false;
		cerr << "ProcessFieldsTD::Process: unknown File-Type" << endl;
	}

	Delete_N_3DArray<FDTD_FLOAT>(field,numLines);

	if (success==false)
	{
		SetEnable(false);
		cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
	}

	return GetNextInterval();
}
