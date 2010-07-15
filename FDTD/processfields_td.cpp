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
#include <iomanip>
#include <sstream>
#include <string>

ProcessFieldsTD::ProcessFieldsTD(Operator* op, Engine* eng) : ProcessFields(op, eng)
{
	pad_length = 8;
}

ProcessFieldsTD::~ProcessFieldsTD()
{
}

void ProcessFieldsTD::DumpNodeInterpol(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling =	Op->GetGridDelta();
#endif

	if (m_DumpType==H_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** H_T = Create_N_3DArray<FDTD_FLOAT>(numLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		unsigned int OpPos[3];
		double delta;
//		cerr << "processing e-fields... " << endl;
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					//in x
					delta  = Op->GetMeshDelta(0,OpPos,true);
					if (delta)
					{
						H_T[0][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(0,OpPos[0],OpPos[1],OpPos[2]) + Eng->GetCurr(0,OpPos[0],OpPos[1]+1,OpPos[2]) + Eng->GetCurr(0,OpPos[0],OpPos[1],OpPos[2]+1) + Eng->GetCurr(0,OpPos[0],OpPos[1]+1,OpPos[2]+1);
						H_T[0][pos[0]][pos[1]][pos[2]] /= (4*delta);
					}
					//in y
					delta  = Op->GetMeshDelta(1,OpPos,true);
					if (delta)
					{
						H_T[1][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(1,OpPos[0],OpPos[1],OpPos[2]) + Eng->GetCurr(1,OpPos[0]+1,OpPos[1],OpPos[2]) + Eng->GetCurr(1,OpPos[0],OpPos[1],OpPos[2]+1) + Eng->GetCurr(1,OpPos[0]+1,OpPos[1],OpPos[2]+1);
						H_T[1][pos[0]][pos[1]][pos[2]] /= (4*delta);//*Op->gridDelta);
					}
					//in z
					delta  = Op->GetMeshDelta(2,OpPos,true);
					if (delta)
					{
						H_T[2][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(2,OpPos[0],OpPos[1],OpPos[2]) + Eng->GetCurr(2,OpPos[0],OpPos[1]+1,OpPos[2]) + Eng->GetCurr(2,OpPos[0]+1,OpPos[1],OpPos[2]) + Eng->GetCurr(2,OpPos[0]+1,OpPos[1]+1,OpPos[2]);
						H_T[2][pos[0]][pos[1]][pos[2]] /= (4*delta);
					}
				}
			}
		}

		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("H-Field"),H_T,discLines,numLines,m_precision,GetInterpolationNameByType(m_DumpMode), m_Mesh_Type, discLines_scaling);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),H_T,numLines,Eng->GetNumberOfTimesteps()*Op->GetTimestep());
		}
		else
			cerr << "ProcessFieldsTD::DumpNodeInterpol: unknown File-Type" << endl;
		Delete_N_3DArray(H_T,numLines);
		H_T = NULL;
	}

	if (m_DumpType==E_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** E_T = Create_N_3DArray<FDTD_FLOAT>(numLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		unsigned int OpPos[3];
		unsigned int OpPosUp[3];
		double delta, deltaUp, deltaRel;
//		cerr << "processing h-fields... " << endl;
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			OpPosUp[0]=start[0]+pos[0]*subSample[0];
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				OpPosUp[1]=start[1]+pos[1]*subSample[1];
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					OpPosUp[2]=start[2]+pos[2]*subSample[2];

					for (int n=0;n<3;++n)
					{
						delta = Op->GetMeshDelta(n,OpPos);
						++OpPosUp[n];
						deltaUp = Op->GetMeshDelta(n,OpPos);
						deltaRel = delta / (delta+deltaUp);
						if (delta*deltaUp)
						{
							E_T[n][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(n,OpPos)*(1-deltaRel)/delta + Eng->GetVolt(n,OpPosUp)/deltaUp*deltaRel;
						}
						--OpPosUp[n];
					}
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("E-Field"),E_T,discLines,numLines,m_precision,GetInterpolationNameByType(m_DumpMode), m_Mesh_Type, discLines_scaling);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),E_T,numLines,(0.5+Eng->GetNumberOfTimesteps())*Op->GetTimestep());
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;
		Delete_N_3DArray(E_T,numLines);
		E_T = NULL;
	}
}

void ProcessFieldsTD::DumpCellInterpol(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = Op->GetGridDelta();
#endif

	if (m_DumpType==E_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** E_T = Create_N_3DArray<FDTD_FLOAT>(numDLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		unsigned int OpPos[3];
		double delta;
//		cerr << "processing e-fields... " << endl;
		for (pos[0]=0;pos[0]<numDLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			for (pos[1]=0;pos[1]<numDLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				for (pos[2]=0;pos[2]<numDLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					//in x
					delta  = Op->GetMeshDelta(0,OpPos); //Op->discLines[0][OpPos[0]+1] - Op->discLines[0][OpPos[0]];
					if (delta)
					{
						E_T[0][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(0,OpPos[0],OpPos[1],OpPos[2]) + Eng->GetVolt(0,OpPos[0],OpPos[1]+1,OpPos[2]) + Eng->GetVolt(0,OpPos[0],OpPos[1],OpPos[2]+1) + Eng->GetVolt(0,OpPos[0],OpPos[1]+1,OpPos[2]+1);
						E_T[0][pos[0]][pos[1]][pos[2]] /= (4*delta);//*Op->gridDelta);
					}
					//in y
					delta  = Op->GetMeshDelta(1,OpPos); //Op->discLines[1][OpPos[1]+1] - Op->discLines[1][OpPos[1]];
					if (delta)
					{
						E_T[1][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(1,OpPos[0],OpPos[1],OpPos[2]) + Eng->GetVolt(1,OpPos[0]+1,OpPos[1],OpPos[2]) + Eng->GetVolt(1,OpPos[0],OpPos[1],OpPos[2]+1) + Eng->GetVolt(1,OpPos[0]+1,OpPos[1],OpPos[2]+1);
						E_T[1][pos[0]][pos[1]][pos[2]] /= (4*delta);//*Op->gridDelta);
					}
					//in z
					delta  = Op->GetMeshDelta(2,OpPos); //Op->discLines[2][OpPos[2]+1] - Op->discLines[2][OpPos[2]];
					if (delta)
					{
						E_T[2][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(2,OpPos[0],OpPos[1],OpPos[2]) + Eng->GetVolt(2,OpPos[0],OpPos[1]+1,OpPos[2]) + Eng->GetVolt(2,OpPos[0]+1,OpPos[1],OpPos[2]) + Eng->GetVolt(2,OpPos[0]+1,OpPos[1]+1,OpPos[2]);
						E_T[2][pos[0]][pos[1]][pos[2]] /= (4*delta);//*Op->gridDelta);
					}
				}
			}
		}

		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("E-Field"),E_T,discDLines,numDLines,m_precision,GetInterpolationNameByType(m_DumpMode), m_Mesh_Type, discLines_scaling);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),E_T,numDLines,Eng->GetNumberOfTimesteps()*Op->GetTimestep());
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;
		Delete_N_3DArray(E_T,numDLines);
		E_T = NULL;
	}

	if (m_DumpType==1)
	{
		//create array
		FDTD_FLOAT**** H_T = Create_N_3DArray<FDTD_FLOAT>(numDLines);
		unsigned int pos[3] = {start[0],start[1],start[2]};
		unsigned int OpPos[3];
		unsigned int OpPosUp[3];
		double delta, deltaUp, deltaRel;
//		cerr << "processing h-fields... " << endl;
		for (pos[0]=0;pos[0]<numDLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			OpPosUp[0]=start[0]+pos[0]*subSample[0];
			for (pos[1]=0;pos[1]<numDLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				OpPosUp[1]=start[1]+pos[1]*subSample[1];
				for (pos[2]=0;pos[2]<numDLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					OpPosUp[2]=start[2]+pos[2]*subSample[2];

					for (int n=0;n<3;++n)
					{
						delta = Op->GetMeshDelta(n,OpPos,true);
						++OpPosUp[n];
						deltaUp = Op->GetMeshDelta(n,OpPos,true);
						deltaRel = delta / (delta+deltaUp);
						if (delta*deltaUp)
						{
							H_T[n][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(n,OpPos)*(1-deltaRel)/delta + Eng->GetCurr(n,OpPosUp)/deltaUp*deltaRel;
						}
						--OpPosUp[n];
					}
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("H-Field"),H_T,discDLines,numDLines,m_precision,GetInterpolationNameByType(m_DumpMode), m_Mesh_Type, discLines_scaling);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),H_T,numDLines,(0.5+Eng->GetNumberOfTimesteps())*Op->GetTimestep());
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;
		Delete_N_3DArray(H_T,numDLines);
		H_T = NULL;
	}
}

void ProcessFieldsTD::DumpNoInterpol(string filename)
{
#ifdef OUTPUT_IN_DRAWINGUNITS
	double discLines_scaling = 1;
#else
	double discLines_scaling = Op->GetGridDelta();
#endif

	unsigned int pos[3];
	unsigned int OpPos[3];
	double delta[3];
	if (m_DumpType==E_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** E_T = Create_N_3DArray<FDTD_FLOAT>(numLines);
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			delta[0]=Op->GetMeshDelta(0,OpPos);//fabs(Op->MainOp->GetIndexDelta(0,OpPos[0]));
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				delta[1]=Op->GetMeshDelta(1,OpPos);//fabs(Op->MainOp->GetIndexDelta(1,OpPos[1]));
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					delta[2]=Op->GetMeshDelta(2,OpPos);//fabs(Op->MainOp->GetIndexDelta(2,OpPos[2]));
					if (delta[0])
						E_T[0][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(0,OpPos[0],OpPos[1],OpPos[2])/delta[0];// /Op->gridDelta;
					if (delta[1])
						E_T[1][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(1,OpPos[0],OpPos[1],OpPos[2])/delta[1];// /Op->gridDelta;
					if (delta[2])
						E_T[2][pos[0]][pos[1]][pos[2]] = Eng->GetVolt(2,OpPos[0],OpPos[1],OpPos[2])/delta[2];// /Op->gridDelta;
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("E-Field"),E_T,discLines,numLines,m_precision,GetInterpolationNameByType(m_DumpMode), m_Mesh_Type, discLines_scaling);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),E_T,numLines,Eng->GetNumberOfTimesteps()*Op->GetTimestep());
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;

		Delete_N_3DArray(E_T,numLines);
		E_T = NULL;
	}

	if (m_DumpType==H_FIELD_DUMP)
	{
		//create array
		FDTD_FLOAT**** H_T = Create_N_3DArray<FDTD_FLOAT>(numLines);
		for (pos[0]=0;pos[0]<numLines[0];++pos[0])
		{
			OpPos[0]=start[0]+pos[0]*subSample[0];
			delta[0]=Op->GetMeshDelta(0,OpPos,true);//fabs(Op->MainOp->GetIndexWidth(0,OpPos[0]));
			for (pos[1]=0;pos[1]<numLines[1];++pos[1])
			{
				OpPos[1]=start[1]+pos[1]*subSample[1];
				delta[1]=Op->GetMeshDelta(1,OpPos,true);//fabs(Op->MainOp->GetIndexWidth(1,OpPos[1]));
				for (pos[2]=0;pos[2]<numLines[2];++pos[2])
				{
					OpPos[2]=start[2]+pos[2]*subSample[2];
					delta[2]=Op->GetMeshDelta(2,OpPos,true);//fabs(Op->MainOp->GetIndexWidth(2,OpPos[2]));
					//in x
					if (delta[0])
						H_T[0][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(0,OpPos[0],OpPos[1],OpPos[2])/delta[0];// /Op->gridDelta;
					if (delta[1])
						H_T[1][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(1,OpPos[0],OpPos[1],OpPos[2])/delta[1];// /Op->gridDelta;
					if (delta[2])
						H_T[2][pos[0]][pos[1]][pos[2]] = Eng->GetCurr(2,OpPos[0],OpPos[1],OpPos[2])/delta[2];// /Op->gridDelta;
				}
			}
		}
		if (m_fileType==VTK_FILETYPE)
		{
			ofstream file(filename.c_str());
			if (file.is_open()==false) { cerr << "ProcessFieldsTD::Process: can't open file '" << filename << "' for writing... abort! " << endl;};
			DumpVectorArray2VTK(file,string("H-Field"),H_T,discLines,numLines,m_precision,GetInterpolationNameByType(m_DumpMode), m_Mesh_Type, discLines_scaling);
			file.close();
		}
		else if (m_fileType==HDF5_FILETYPE)
		{
			stringstream ss;
			ss << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();
			DumpVectorArray2HDF5(filename.c_str(),string( ss.str() ),H_T,numLines,(0.5+Eng->GetNumberOfTimesteps())*Op->GetTimestep());
		}
		else
			cerr << "ProcessFieldsTD::DumpCellInterpol: unknown File-Type" << endl;

		Delete_N_3DArray(H_T,numLines);
		H_T = NULL;
	}
}

int ProcessFieldsTD::Process()
{
	if (Enabled==false) return -1;
	if (filePattern.empty()) return -1;
	if (CheckTimestep()==false) return GetNextInterval();
	stringstream ss;
	ss << filePattern << std::setw( pad_length ) << std::setfill( '0' ) << Eng->GetNumberOfTimesteps();

	if (m_fileType==VTK_FILETYPE)
	{
		ss << ".vtk";
		if (m_DumpMode==NO_INTERPOLATION)
			DumpNoInterpol(ss.str());
		if (m_DumpMode==NODE_INTERPOLATE)
			DumpNodeInterpol(ss.str());
		if (m_DumpMode==CELL_INTERPOLATE)
			DumpCellInterpol(ss.str());
	}
	else if (m_fileType==HDF5_FILETYPE)
	{
		if (m_DumpMode==NO_INTERPOLATION)
			DumpNoInterpol(m_filename);
		if (m_DumpMode==NODE_INTERPOLATE)
			DumpNodeInterpol(m_filename);
		if (m_DumpMode==CELL_INTERPOLATE)
			DumpCellInterpol(m_filename);
	}
	return GetNextInterval();
}
