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
public:
    Operator();
	virtual ~Operator();

	virtual void SetGeometryCSX(ContinuousStructure* geo);

	virtual int CalcECOperator() {};

	virtual void ApplyElectricBC(bool* dirs) {}; //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs) {};

	double GetTimestep() {return dT;};

	virtual void Reset();

protected:
	virtual void Init();
	ContinuousStructure* CSX;
	double gridDelta;

	double* discLines[3];
	unsigned int numLines[3];

	//EC operator
	FDTD_FLOAT* vv[3]; //calc new voltage from old voltage
	FDTD_FLOAT* vi[3]; //calc new voltage from old current
	FDTD_FLOAT* ii[3]; //calc new current from old current
	FDTD_FLOAT* iv[3]; //calc new current from old voltage

	//E-Field Excitation
	//!	  Calc the electric field excitation.
	virtual bool CalcEFieldExcitation() {};
	unsigned int E_Ex_Count;
	unsigned int* E_Ex_index;
	FDTD_FLOAT* E_Ex_amp[3]; //represented as edge-voltages!!
	FDTD_FLOAT* E_Ex_delay;

	//Calc timestep only internal use
	virtual double CalcTimestep() {};
	double dT; //FDTD timestep!
};

#endif // OPERATOR_H
