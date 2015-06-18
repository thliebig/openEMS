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
#include "tools/sar_calculation.h"

#include "CSPropMaterial.h"

using namespace std;

ProcessFieldsSAR::ProcessFieldsSAR(Engine_Interface_Base* eng_if) : ProcessFieldsFD(eng_if)
{
	m_UseCellKappa = true;
	m_SAR_method = "Simple";
}

ProcessFieldsSAR::~ProcessFieldsSAR()
{
	for (size_t n = 0; n<m_E_FD_Fields.size(); ++n)
		Delete_N_3DArray(m_E_FD_Fields.at(n),numLines);
	m_E_FD_Fields.clear();

	for (size_t n = 0; n<m_J_FD_Fields.size(); ++n)
		Delete_N_3DArray(m_J_FD_Fields.at(n),numLines);
	m_J_FD_Fields.clear();
}

void ProcessFieldsSAR::SetDumpType(DumpType type)
{
	if (type==SAR_RAW_DATA)
		m_UseCellKappa = true;
	ProcessFieldsFD::SetDumpType(type);
}

bool ProcessFieldsSAR::NeedConductivity() const
{
	return !m_UseCellKappa;
}

void ProcessFieldsSAR::SetSubSampling(unsigned int subSampleRate, int dir)
{
	UNUSED(subSampleRate);UNUSED(dir);
	cerr << "ProcessFieldsSAR::SetSubSampling: Warning: Defining a sub-sampling for SAR calculation is not allowed!!! Skipped!" << endl;
}

void ProcessFieldsSAR::SetOptResolution(double optRes, int dir)
{
	UNUSED(optRes);UNUSED(dir);
	cerr << "ProcessFieldsSAR::SetOptResolution: Warning: Defining a sub-sampling for SAR calculation is not allowed!!! Skipped!" << endl;
}

void ProcessFieldsSAR::InitProcess()
{
	if ((m_DumpType!=SAR_LOCAL_DUMP) && (m_DumpType!=SAR_1G_DUMP) && (m_DumpType!=SAR_10G_DUMP) && (m_DumpType!=SAR_RAW_DATA))
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

	if ((m_DumpType==SAR_RAW_DATA) && (m_fileType!=HDF5_FILETYPE))
	{
		Enabled=false;
		cerr << "ProcessFieldsSAR::InitProcess(): Error, wrong file type for dumping raw SAR data! skipping" << endl;
		return;

	}

	ProcessFieldsFD::InitProcess();

	if (Enabled==false) return;

	//create data structures...
	for (size_t n = 0; n<m_FD_Samples.size(); ++n)
	{
		m_E_FD_Fields.push_back(Create_N_3DArray<std::complex<float> >(numLines));
		if (!m_UseCellKappa)
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

	//save dump type
	DumpType save_dump_type = m_DumpType;

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
	if (!m_UseCellKappa)
	{
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
	}

	//reset dump type
	m_DumpType = save_dump_type;

	++m_FD_SampleCount;
	return GetNextInterval();
}

void ProcessFieldsSAR::DumpFDData()
{
	if (Enabled==false) return;
	unsigned int pos[3];
	unsigned int orig_pos[3];
	float*** SAR = Create3DArray<float>(numLines);
	double coord[3];
	ContinuousStructure* CSX = Op->GetGeometryCSX();
	CSProperties* prop = NULL;
	CSPropMaterial* matProp = NULL;

	double power;

	float*** cell_volume = Create3DArray<float>(numLines);
	float*** cell_density = Create3DArray<float>(numLines);
	float*** cell_kappa = NULL;
	if (m_UseCellKappa)
		cell_kappa = Create3DArray<float>(numLines);

	bool found_UnIsotropic=false;

	// calculate volumes and masses for all cells
	for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
	{
		orig_pos[0] = posLines[0][pos[0]];
		for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
		{
			orig_pos[1] = posLines[1][pos[1]];
			vector<CSPrimitives*> vPrims = Op->GetPrimitivesBoundBox(orig_pos[0], orig_pos[1], -1, CSProperties::MATERIAL);
			for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
			{
				orig_pos[2] = posLines[2][pos[2]];

				cell_volume[pos[0]][pos[1]][pos[2]] = Op->GetCellVolume(orig_pos);
				cell_density[pos[0]][pos[1]][pos[2]] = 0.0;

				Op->GetCellCenterMaterialAvgCoord(orig_pos, coord);
				prop = CSX->GetPropertyByCoordPriority(coord, vPrims);
//				prop = CSX->GetPropertyByCoordPriority(coord,CSProperties::MATERIAL);
				if (prop!=0)
				{
					matProp = dynamic_cast<CSPropMaterial*>(prop);
					if (matProp)
					{
						found_UnIsotropic |= !matProp->GetIsotropy();
						cell_density[pos[0]][pos[1]][pos[2]] = matProp->GetDensityWeighted(coord);
						if (m_UseCellKappa)
							cell_kappa[pos[0]][pos[1]][pos[2]] = matProp->GetKappaWeighted(0,coord);
					}
				}
			}
		}
	}
	if (found_UnIsotropic)
		cerr << "ProcessFieldsSAR::DumpFDData(): Warning, found unisotropic material in SAR calculation... this is unsupported!" << endl;

	float* cellWidth[3];
	for (int n=0;n<3;++n)
	{
		cellWidth[n]=new float[numLines[n]];
		for (unsigned int i=0;i<numLines[n];++i)
			cellWidth[n][i]=Op->GetDiscDelta(n,posLines[n][i])*Op->GetGridDelta();
	}

	if (m_DumpType == SAR_RAW_DATA)
	{
		if (m_fileType!=HDF5_FILETYPE)
		{
			cerr << "ProcessFieldsSAR::DumpFDData(): Error, wrong file type, this should not happen!!! skipped" << endl;
			return;
		}

		size_t datasize[]={numLines[0],numLines[1],numLines[2]};
		for (size_t n = 0; n<m_FD_Samples.size(); ++n)
		{
			stringstream ss;
			ss << "f" << n;
			if (m_HDF5_Dump_File->WriteVectorField(ss.str(), m_E_FD_Fields.at(n), datasize)==false)
				cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
		}

		m_HDF5_Dump_File->SetCurrentGroup("/CellData");
		if (m_UseCellKappa==false)
			cerr <<  "ProcessFieldsSAR::DumpFDData: Error, cell conductivity data not available, this should not happen... skipping! " << endl;
		else if (m_HDF5_Dump_File->WriteScalarField("Conductivity", cell_kappa, datasize)==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
		if (m_HDF5_Dump_File->WriteScalarField("Density", cell_density, datasize)==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
		if (m_HDF5_Dump_File->WriteScalarField("Volume", cell_volume, datasize)==false)
			cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
	}
	else
	{
		SAR_Calculation SAR_Calc;
		SAR_Calc.SetAveragingMethod(m_SAR_method, g_settings.GetVerboseLevel()==0);
		SAR_Calc.SetDebugLevel(g_settings.GetVerboseLevel());
		SAR_Calc.SetNumLines(numLines);
		if (m_DumpType == SAR_LOCAL_DUMP)
			SAR_Calc.SetAveragingMass(0);
		else if (m_DumpType == SAR_1G_DUMP)
			SAR_Calc.SetAveragingMass(1e-3);
		else if (m_DumpType == SAR_10G_DUMP)
			SAR_Calc.SetAveragingMass(10e-3);
		else
		{
			cerr << "ProcessFieldsSAR::DumpFDData: unknown SAR dump type...!" << endl;
		}
		SAR_Calc.SetCellDensities(cell_density);
		SAR_Calc.SetCellWidth(cellWidth);
		SAR_Calc.SetCellVolumes(cell_volume);
		SAR_Calc.SetCellCondictivity(cell_kappa); // cell_kappa will be NULL if m_UseCellKappa is false

		for (size_t n = 0; n<m_FD_Samples.size(); ++n)
		{
			SAR_Calc.SetEField(m_E_FD_Fields.at(n));
			if (!m_UseCellKappa)
				SAR_Calc.SetJField(m_J_FD_Fields.at(n));
			power = SAR_Calc.CalcSARPower();
			SAR_Calc.CalcSAR(SAR);

			if (m_fileType==VTK_FILETYPE)
			{
				stringstream ss;
				ss << m_filename << fixed << "_f=" << m_FD_Samples.at(n);

				m_Vtk_Dump_File->SetFilename(ss.str());
				m_Vtk_Dump_File->ClearAllFields();
				m_Vtk_Dump_File->AddScalarField(GetFieldNameByType(m_DumpType),SAR);
				if (m_Vtk_Dump_File->Write()==false)
					cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
			}
			else if (m_fileType==HDF5_FILETYPE)
			{
				stringstream ss;
				ss << "f" << n;
				size_t datasize[]={numLines[0],numLines[1],numLines[2]};
				if (m_HDF5_Dump_File->WriteScalarField(ss.str(), SAR, datasize)==false)
					cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
				float freq[1] = {(float)m_FD_Samples.at(n)};
				if (m_HDF5_Dump_File->WriteAtrribute("/FieldData/FD/"+ss.str(),"frequency",freq,1)==false)
					cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
				float pow[1] = {(float)power};
				if (m_HDF5_Dump_File->WriteAtrribute("/FieldData/FD/"+ss.str(),"power",pow,1)==false)
					cerr << "ProcessFieldsSAR::DumpFDData: can't dump to file...! " << endl;
			}
			else
				cerr << "ProcessFieldsSAR::DumpFDData: unknown File-Type" << endl;
		}
	}
	for (int n=0;n<3;++n)
		delete[] cellWidth[n];
	Delete3DArray(cell_volume,numLines);
	Delete3DArray(cell_density,numLines);
	Delete3DArray(cell_kappa,numLines);
	Delete3DArray(SAR,numLines);
}
