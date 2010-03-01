#include "operator.h"

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

	E_Ex_index = NULL;
	E_Ex_delay = NULL;
	for (int n=0;n<3;++n)
	{
		discLines[n]=NULL;
		vv[n]=NULL;
		vi[n]=NULL;
		iv[n]=NULL;
		ii[n]=NULL;
		E_Ex_amp[n]=NULL;
	}
}

void Operator::Reset()
{
	delete[] E_Ex_index;
	delete[] E_Ex_delay;
	for (int n=0;n<3;++n)
	{
		delete[] vv[n];
		delete[] vi[n];
		delete[] iv[n];
		delete[] ii[n];
		delete[] E_Ex_amp[n];
	}
	Operator::Init();
}

void Operator::SetGeometryCSX(ContinuousStructure* geo)
{
	if (geo==NULL) return;

	Reset();
	CSX = geo;
}
