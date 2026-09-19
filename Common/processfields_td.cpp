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
#include "async_dumps.h"
#include <iomanip>
#include <sstream>
#include <string>

using namespace std;

ProcessFieldsTD::ProcessFieldsTD(Engine_Interface_Base* eng_if) : ProcessFields(eng_if)
{
	pad_length = 8;
	m_AsyncFailed = false;
	m_AsyncUsed = false;
	m_Snapshots = false;
}

ProcessFieldsTD::~ProcessFieldsTD()
{
	FinishAsync();
}

void ProcessFieldsTD::FinishAsync()
{
	if (!m_AsyncUsed)
		return;
	AsyncDumps::Get().Wait(this);
	m_AsyncUsed = false;
	if (m_AsyncFailed)
	{
		SetEnable(false);
		cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
		m_AsyncFailed = false;
	}
}

// on the background writer: calculate the field from the snapshot \a src (if \a field is NULL) and write it
bool ProcessFieldsTD::WriteHDF5(unsigned int ts, float time, ArrayLib::ArrayNIJK<float>* field, const float* src)
{
	ArrayLib::ArrayNIJK<float> calc_field;
	if (!field)
	{
		calc_field.Init("TD_field", numLines);
		if (!CalcField(calc_field, src))
			return false;
		field = &calc_field;
	}
	stringstream ss;
	ss << std::setw( pad_length ) << std::setfill( '0' ) << ts;
	bool success = m_HDF5_Dump_File->WriteVectorField<float>(ss.str(), *field, g_settings.GetLegacyHDF5Dumps());
	float time_attr[1] = {time};
	success &= m_HDF5_Dump_File->WriteAttribute("/FieldData/TD/"+ss.str(),"time",time_attr,1);
	return success;
}

void ProcessFieldsTD::InitProcess()
{
	m_Snapshots = false;
	if (Enabled==false) return;

	ProcessFields::InitProcess();

	// before the first snapshot, e.g. to evaluate the dumped nodes on the device
	if (m_fileType==HDF5_FILETYPE)
	{
		GetGather();
		m_Snapshots = m_Gather && m_Eng_Interface->PrepareSnapshotGather(m_Gather);
	}

	if (m_Vtk_Dump_File)
		m_Vtk_Dump_File->SetHeader(string("openEMS TD Field Dump -- Interpolation: ")+m_Eng_Interface->GetInterpolationTypeString());

	if (m_HDF5_Dump_File)
	{
		m_HDF5_Dump_File->SetCurrentGroup("/FieldData/TD");
		// a dataset every few timesteps: keep the file open, see PostProcess()
		m_HDF5_Dump_File->SetKeepOpen(true);
	}
}

void ProcessFieldsTD::PostProcess()
{
	FinishAsync();
	ProcessFields::PostProcess();
	if (m_HDF5_Dump_File)
		m_HDF5_Dump_File->Close();
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (m_AsyncFailed)
	{
		FinishAsync();   // reports and disables
		return -1;
	}
	if (CheckTimestep()==false) return GetNextInterval();

	string filename = m_filename;

	// HDF5: written in the background, see AsyncDumps; with a snapshot of the fields also calculated there
	if (m_fileType==HDF5_FILETYPE)
	{
		const unsigned int ts = m_Eng_Interface->GetNumberOfTimesteps();
		const float time = (float)m_Eng_Interface->GetTime(m_dualTime);
		int slot = -1;
		const float *volt = NULL, *curr = NULL;
		const float* src = NULL;
		if (m_Snapshots && AsyncDumps::Get().Snapshot(m_Eng_Interface, ts, slot, volt, curr))
			src = (m_DumpType==H_FIELD_DUMP) ? curr : volt;
		ArrayLib::ArrayNIJK<float>* field = NULL;
		if (!src)
		{
			field = new ArrayLib::ArrayNIJK<float>("TD_field", numLines);
			if (!CalcField(*field))
			{
				delete field;
				SetEnable(false);
				cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
				return GetNextInterval();
			}
		}
		m_AsyncUsed = true;
		AsyncDumps::Get().Push(this, slot, [this, ts, time, field, src, slot]()
		{
			if (src)
				m_Eng_Interface->WaitFieldSnapshot(slot);   // the copy may still run on the device
			if (!WriteHDF5(ts, time, field, src))
				m_AsyncFailed = true;
			delete field;
		});
		return GetNextInterval();
	}

	ArrayLib::ArrayNIJK<float> field("TD_field", numLines);
	bool success = CalcField(field);

	if (m_fileType==VTK_FILETYPE)
	{
		m_Vtk_Dump_File->SetTimestep(m_Eng_Interface->GetNumberOfTimesteps());
		m_Vtk_Dump_File->ClearAllFields();
		m_Vtk_Dump_File->AddVectorField(GetFieldNameByType(m_DumpType),field);
		success &= m_Vtk_Dump_File->Write();
	}
	else
	{
		success = false;
		cerr << "ProcessFieldsTD::Process: unknown File-Type" << endl;
	}

	if (success==false)
	{
		SetEnable(false);
		cerr << "ProcessFieldsTD::Process: can't dump to file... disabled! " << endl;
	}

	return GetNextInterval();
}
