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
}

void ProcessFieldsSAR::InitProcess()
{
	if (m_DumpType!=SAR_LOCAL_DUMP)
	{
		Enabled=false;
		cerr << "ProcessFieldsSAR::InitProcess(): Error, wrong dump type... this should not happen... skipping!" << endl;
		return;
	}
	if (m_Eng_Interface->GetInterpolationType()!=Engine_Interface_Base::NODE_INTERPOLATE)
	{
		cerr << "ProcessFieldsSAR::InitProcess(): Warning, interpolation type is not supported, resetting to NODE!" << endl;
		SetDumpMode2Node();
	}
	ProcessFieldsFD::InitProcess();
}

double ProcessFieldsSAR::GetKappaDensityRatio(const unsigned int* pos)
{

	double coord[3] = {discLines[0][pos[0]],discLines[1][pos[1]],discLines[2][pos[2]]};
	unsigned int OpPos[] = {posLines[0][pos[0]],posLines[1][pos[1]],posLines[2][pos[2]]};
	ContinuousStructure* CSX = Op->GetGeometryCSX();
	CSProperties* prop = CSX->GetPropertyByCoordPriority(coord,CSProperties::MATERIAL);

	if (prop==0)
		return 0.0;
	CSPropMaterial* matProp = dynamic_cast<CSPropMaterial*>(prop);
	double density = matProp->GetDensityWeighted(coord);

	if (density==0)
		return 0.0;

	double kappa = 0;
	double max_kappa = 0;

	for (int n=0;n<3;++n)
	{
		kappa = Op->GetDiscMaterial(1,n,OpPos);
		if (kappa>max_kappa)
			max_kappa = kappa;
		if (OpPos[n]>0)
		{
			--OpPos[n];
			kappa = Op->GetDiscMaterial(1,n,OpPos);
			if (kappa>max_kappa)
				max_kappa = kappa;
			++OpPos[n];
		}
	}

	return max_kappa/density;
}

void ProcessFieldsSAR::DumpFDData()
{
	unsigned int pos[3];
	FDTD_FLOAT*** SAR = Create3DArray<float>(numLines);
	std::complex<float>**** field_fd = NULL;
	double kdRatio = 0;

	for (size_t n = 0; n<m_FD_Samples.size(); ++n)
	{
		field_fd = m_FD_Fields.at(n);
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
		{
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
			{
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					kdRatio = GetKappaDensityRatio(pos);
					SAR[pos[0]][pos[1]][pos[2]] = kdRatio*pow(abs(field_fd[0][pos[0]][pos[1]][pos[2]]) , 2);
					SAR[pos[0]][pos[1]][pos[2]] = kdRatio*pow(abs(field_fd[1][pos[0]][pos[1]][pos[2]]) , 2);
					SAR[pos[0]][pos[1]][pos[2]] = kdRatio*pow(abs(field_fd[2][pos[0]][pos[1]][pos[2]]) , 2);
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
			float freq[1]={m_FD_Samples.at(n)};
			if (m_HDF5_Dump_File->WriteAtrribute("/FieldData/FD/"+ss.str(),"frequency",freq,1)==false)
				cerr << "ProcessFieldsSAR::Process: can't dump to file...! " << endl;
		}
		else
			cerr << "ProcessFieldsSAR::Process: unknown File-Type" << endl;
	}

	Delete3DArray(SAR,numLines);
}
