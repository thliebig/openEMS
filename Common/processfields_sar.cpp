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

#include "processfields_sar.h"
#include "operator_base.h"
#include "tools/vtk_file_writer.h"
#include "tools/hdf5_file_writer.h"

ProcessFieldsSAR::ProcessFieldsSAR(Engine_Interface_Base* eng_if) : ProcessFieldsFD(eng_if)
{
}

ProcessFieldsSAR::~ProcessFieldsSAR()
{
	for (size_t n = 0; n<m_FD_Fields.size(); ++n)
	{
		Delete_N_3DArray(m_E_FD_Fields.at(n),numLines);
		Delete_N_3DArray(m_J_FD_Fields.at(n),numLines);
	}
	m_E_FD_Fields.clear();
	m_J_FD_Fields.clear();
}

void ProcessFieldsSAR::SetSubSampling(unsigned int subSampleRate, int dir)
{
	cerr << "ProcessFieldsSAR::SetSubSampling: Warning: Defining a sub-sampling for SAR calculation is not recommended!!! Expect false results!" << endl;
	ProcessFieldsFD::SetSubSampling(subSampleRate,dir);
}

void ProcessFieldsSAR::SetOptResolution(double optRes, int dir)
{
	cerr << "ProcessFieldsSAR::SetOptResolution: Warning: Defining a sub-sampling for SAR calculation is not recommended!!! Expect false results!" << endl;
	ProcessFieldsFD::SetOptResolution(optRes,dir);
}

void ProcessFieldsSAR::InitProcess()
{
	if (m_DumpType!=SAR_LOCAL_DUMP)
	{
		Enabled=false;
		cerr << "ProcessFieldsSAR::InitProcess(): Error, wrong dump type... this should not happen... skipping!" << endl;
		return;
	}
	if (m_Eng_Interface->GetInterpolationType()!=Engine_Interface_Base::CELL_INTERPOLATE)
	{
		cerr << "ProcessFieldsSAR::InitProcess(): Warning, interpolation type is not supported, resetting to CELL!" << endl;
		SetDumpMode2Cell();
	}
	ProcessFieldsFD::InitProcess();

	if (Enabled==false) return;

	//create data structures...
	for (size_t n = 0; n<m_FD_Samples.size(); ++n)
	{
		m_E_FD_Fields.push_back(Create_N_3DArray<std::complex<float> >(numLines));
		m_J_FD_Fields.push_back(Create_N_3DArray<std::complex<float> >(numLines));
	}
}

int ProcessFieldsSAR::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();

	if ((m_FD_Interval==0) || (m_Eng_Interface->GetNumberOfTimesteps()%m_FD_Interval!=0))
		return GetNextInterval();

	std::complex<float>**** field_fd = NULL;
	unsigned int pos[3];
	double T;
	FDTD_FLOAT**** field_td=NULL;

	// calc E-field
	m_DumpType = E_FIELD_DUMP;
	field_td = CalcField();
	T = m_Eng_Interface->GetTime(m_dualTime);
	for (size_t n = 0; n<m_FD_Samples.size(); ++n)
	{
		std::complex<float> exp_jwt_2_dt = std::exp( (std::complex<float>)(-2.0 * _I * M_PI * m_FD_Samples.at(n) * T) );
		exp_jwt_2_dt *= 2; // *2 for single-sided spectrum
		exp_jwt_2_dt *= Op->GetTimestep() * m_FD_Interval; // multiply with timestep-interval
		field_fd = m_E_FD_Fields.at(n);
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					field_fd[0][pos[0]][pos[1]][pos[2]] += field_td[0][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
					field_fd[1][pos[0]][pos[1]][pos[2]] += field_td[1][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
					field_fd[2][pos[0]][pos[1]][pos[2]] += field_td[2][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
				}
			}
		}
	}
	Delete_N_3DArray<FDTD_FLOAT>(field_td,numLines);

	// calc J-field
	m_DumpType = J_FIELD_DUMP;
	field_td = CalcField();
	T = m_Eng_Interface->GetTime(m_dualTime);
	for (size_t n = 0; n<m_FD_Samples.size(); ++n)
	{
		std::complex<float> exp_jwt_2_dt = std::exp( (std::complex<float>)(-2.0 * _I * M_PI * m_FD_Samples.at(n) * T) );
		exp_jwt_2_dt *= 2; // *2 for single-sided spectrum
		exp_jwt_2_dt *= Op->GetTimestep() * m_FD_Interval; // multiply with timestep-interval
		field_fd = m_J_FD_Fields.at(n);
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					field_fd[0][pos[0]][pos[1]][pos[2]] += field_td[0][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
					field_fd[1][pos[0]][pos[1]][pos[2]] += field_td[1][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
					field_fd[2][pos[0]][pos[1]][pos[2]] += field_td[2][pos[0]][pos[1]][pos[2]] * exp_jwt_2_dt;
				}
			}
		}
	}
	Delete_N_3DArray<FDTD_FLOAT>(field_td,numLines);

	//reset dump type
	m_DumpType = SAR_LOCAL_DUMP;

	++m_FD_SampleCount;
	return GetNextInterval();
}

void ProcessFieldsSAR::DumpFDData()
{
	if (Enabled==false) return;
	unsigned int pos[3];
	unsigned int orig_pos[3];
	FDTD_FLOAT*** SAR = Create3DArray<float>(numLines);
	std::complex<float>**** E_field_fd = NULL;
	std::complex<float>**** J_field_fd = NULL;
	double coord[3];
	double density;
	ContinuousStructure* CSX = Op->GetGeometryCSX();
	CSProperties* prop = NULL;
	CSPropMaterial* matProp = NULL;

	double power;
	double l_pow;

	for (size_t n = 0; n<m_FD_Samples.size(); ++n)
	{
		E_field_fd = m_E_FD_Fields.at(n);
		J_field_fd = m_J_FD_Fields.at(n);
		power = 0;
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			orig_pos[0] = posLines[0][pos[0]];
			coord[0] = Op->GetDiscLine(0,orig_pos[0],true);
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				orig_pos[1] = posLines[1][pos[1]];
				coord[1] = Op->GetDiscLine(1,orig_pos[1],true);
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					orig_pos[2] = posLines[2][pos[2]];
					coord[2] = Op->GetDiscLine(2,orig_pos[2],true);
					l_pow  = abs(E_field_fd[0][pos[0]][pos[1]][pos[2]]) * abs(J_field_fd[0][pos[0]][pos[1]][pos[2]]);
					l_pow += abs(E_field_fd[1][pos[0]][pos[1]][pos[2]]) * abs(J_field_fd[1][pos[0]][pos[1]][pos[2]]);
					l_pow += abs(E_field_fd[2][pos[0]][pos[1]][pos[2]]) * abs(J_field_fd[2][pos[0]][pos[1]][pos[2]]);

					power += 0.5*l_pow*Op->GetCellVolume(orig_pos);

					prop = CSX->GetPropertyByCoordPriority(coord,CSProperties::MATERIAL);
					SAR[pos[0]][pos[1]][pos[2]] = 0.0;
					density=0.0;
					if (prop!=0)
					{
						matProp = dynamic_cast<CSPropMaterial*>(prop);
						density = matProp->GetDensityWeighted(coord);
						if (density>0)
						{
							SAR[pos[0]][pos[1]][pos[2]] = l_pow*0.5/density;
						}
					}
				}
			}
		}

		if (m_fileType==VTK_FILETYPE)
		{
			stringstream ss;
			ss << m_filename << fixed << "_f=" << m_FD_Samples.at(n);

			m_Vtk_Dump_File->SetFilename(ss.str());
			m_Vtk_Dump_File->ClearAllFields();
			m_Vtk_Dump_File->AddScalarField(GetFieldNameByType(m_DumpType),SAR);
			if (m_Vtk_Dump_File->Write()==false)
				cerr << "ProcessFieldsSAR::Process: can't dump to file...! " << endl;
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << "f" << n;
			size_t datasize[]={numLines[0],numLines[1],numLines[2]};
			if (m_HDF5_Dump_File->WriteScalarField(ss.str(), SAR, datasize)==false)
				cerr << "ProcessFieldsSAR::Process: can't dump to file...! " << endl;
			float freq[1] = {(float)m_FD_Samples.at(n)};
			if (m_HDF5_Dump_File->WriteAtrribute("/FieldData/FD/"+ss.str(),"frequency",freq,1)==false)
				cerr << "ProcessFieldsSAR::Process: can't dump to file...! " << endl;
			float pow[1] = {(float)power};
			if (m_HDF5_Dump_File->WriteAtrribute("/FieldData/FD/"+ss.str(),"power",pow,1)==false)
				cerr << "ProcessFieldsSAR::Process: can't dump to file...! " << endl;
		}
		else
			cerr << "ProcessFieldsSAR::Process: unknown File-Type" << endl;
	}

	Delete3DArray(SAR,numLines);
}
