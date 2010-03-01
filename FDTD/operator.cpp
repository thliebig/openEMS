#include <fstream>
#include "operator.h"
#include "tools/array_ops.h"

Operator::Operator()
{
	Operator::Init();
}

Operator::~Operator()
{
	Operator::Reset();
}

void Operator::Init()
{
	CSX = NULL;

	ExciteSignal = NULL;
	E_Ex_delay = NULL;
	vv=NULL;
	vi=NULL;
	iv=NULL;
	ii=NULL;
	for (int n=0;n<3;++n)
	{
		discLines[n]=NULL;
		E_Ex_amp[n]=NULL;
		E_Ex_index[n]=NULL;
	}
}

void Operator::Reset()
{
	delete[] ExciteSignal;
	delete[] E_Ex_delay;
	Delete_N_3DArray(vv,numLines);
	Delete_N_3DArray(vi,numLines);
	Delete_N_3DArray(iv,numLines);
	Delete_N_3DArray(ii,numLines);
	for (int n=0;n<3;++n)
	{
		delete[] discLines[n];
		delete[] E_Ex_amp[n];
		delete[] E_Ex_index[n];
	}
	Operator::Init();
}

unsigned int Operator::GetNyquistNum(double fmax)
{
	if (dT==0) return 1;
	double T0 = 1/fmax;
	return floor(T0/2/dT);
}

bool Operator::SnapToMesh(double* dcoord, unsigned int* uicoord)
{
	bool ok=true;
	for (int n=0;n<3;++n)
	{
		if (dcoord[n]<discLines[n][0]) {ok=false;uicoord[n]=0;};
		if (dcoord[n]>discLines[n][numLines[n]-1]) {ok=false;uicoord[n]=numLines[n]-1;};
		for (unsigned int i=0;i<numLines[n]-1;++i)
		{
			if (dcoord[n]<=discLines[n][i])
			{
				if (fabs(dcoord[n]-discLines[n][i])<(fabs(dcoord[n]-discLines[n][i+1])))
					uicoord[n]=i;
				else
					uicoord[n]=i+1;
				i = numLines[n];
			}
		}
	}
//	cerr << "Operator::SnapToMesh: " << discLines[0][uicoord[0]] << " " << discLines[1][uicoord[1]] << " " << discLines[2][uicoord[2]] << endl;
//	cerr << "Operator::SnapToMesh: " << uicoord[0] << " " << uicoord[1] << " " << uicoord[2] << endl;
	return ok;
}

void Operator::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return;

	Reset();
	CSX = geo;
}

double Operator::GetNumberCells()
{
	if (numLines)
		return (numLines[0]-1)*(numLines[1]-1)*(numLines[2]-1);
	return 0;
}

void Operator::ShowSize()
{
	unsigned int OpSize = 12*numLines[0]*numLines[1]*numLines[2]*sizeof(FDTD_FLOAT);
	unsigned int FieldSize = 6*numLines[0]*numLines[1]*numLines[2]*sizeof(FDTD_FLOAT);
	double MBdiff = 1024*1024;

	cout << "FDTD Operator Size:" << endl;
	cout << "Size of Operator in Byte  : " << OpSize << " Byte (" << (double)OpSize/MBdiff << " MB) " << endl;
	cout << "Size of Field-Data in Byte: " << FieldSize << " Byte (" << (double)FieldSize/MBdiff << " MB) " << endl;
}

void Operator::CalcGaussianPulsExcitation(double f0, double fc)
{
	if (dT==0) return;

	ExciteLength = (unsigned int)(2.0 * 9.0/(2.0*PI*fc) / dT);
	cerr << "Operator::CalcGaussianPulsExcitation: Length of the excite signal: " << ExciteLength << " timesteps" << endl;
	delete[] ExciteSignal;
	ExciteSignal = new FDTD_FLOAT[ExciteLength+1];
	ExciteSignal[0]=0.0;
	for (unsigned int n=1;n<ExciteLength+1;++n)
	{
		ExciteSignal[n] = cos(2.0*PI*f0*(n*dT-9.0/(2.0*PI*fc)))*exp(-1*pow(2.0*PI*fc*n*dT/3.0-3,2));
//		cerr << ExciteSignal[n] << endl;
	}
}

void Operator::DumpOperator2File(string filename)
{
	ofstream file(filename.c_str(),ios_base::out);
//	file.open;
	if (file.is_open()==false)
	{
		cerr << "Operator::DumpOperator2File: Can't open file: " << filename << endl;
		return;
	}
	file << "########### Operator vv ###########" << endl;
	file << "ix\tiy\tiz\tvv_x\tvv_y\tvv_z" << endl;
	Dump_N_3DArray2File(file,vv,numLines);

	file << "########### Operator vi ###########" << endl;
	file << "ix\tiy\tiz\tvi_x\tvi_y\tvi_z" << endl;
	Dump_N_3DArray2File(file,vi,numLines);

	file << "########### Operator iv ###########" << endl;
	file << "ix\tiy\tiz\tiv_x\tiv_y\tiv_z" << endl;
	Dump_N_3DArray2File(file,iv,numLines);

	file << "########### Operator ii ###########" << endl;
	file << "ix\tiy\tiz\tii_x\tii_y\tii_z" << endl;
	Dump_N_3DArray2File(file,ii,numLines);

	file.close();
}
