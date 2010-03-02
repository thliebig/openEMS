#ifndef OPERATOR_H
#define OPERATOR_H

#include "ContinuousStructure.h"
#include "tools/AdrOp.h"
#include "tools/constants.h"

#define FDTD_FLOAT float

//! Abstract base-class for the FDTD-operator
class Operator
{
	friend class Engine;
	friend class ProcessFields;
	friend class ProcessFieldsTD;
public:
    Operator();
	virtual ~Operator();

	virtual void SetGeometryCSX(ContinuousStructure* geo);

	virtual int CalcECOperator() {return -1;};

	virtual void CalcGaussianPulsExcitation(double f0, double fc);

	virtual void ApplyElectricBC(bool* dirs) {}; //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs) {};

	double GetTimestep() {return dT;};
	unsigned int GetNyquistNum(double fmax);
	double GetNumberCells();

	void ShowSize();

	void DumpOperator2File(string filename);

	virtual void Reset();

	bool SnapToMesh(double* coord, unsigned int* uicoord, bool lower=false);

protected:
	virtual void Init();
	ContinuousStructure* CSX;
	double gridDelta;

	double* discLines[3];
	unsigned int numLines[3];

	//EC operator
	FDTD_FLOAT**** vv; //calc new voltage from old voltage
	FDTD_FLOAT**** vi; //calc new voltage from old current
	FDTD_FLOAT**** ii; //calc new current from old current
	FDTD_FLOAT**** iv; //calc new current from old voltage

	//Excitation time-signal
	unsigned int ExciteLength;
	FDTD_FLOAT* ExciteSignal;

	//E-Field Excitation
	//!	  Calc the electric field excitation.
	virtual bool CalcEFieldExcitation() {return false;};
	unsigned int E_Ex_Count;
	unsigned int* E_Ex_index[3];
	FDTD_FLOAT* E_Ex_amp[3]; //represented as edge-voltages!!
	unsigned int* E_Ex_delay;

	//Calc timestep only internal use
	virtual double CalcTimestep() {return 0;};
	double dT; //FDTD timestep!
};

#endif // OPERATOR_H
