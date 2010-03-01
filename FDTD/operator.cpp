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
