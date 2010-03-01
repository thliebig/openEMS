#ifndef CARTOPERATOR_H
#define CARTOPERATOR_H

#include "ContinuousStructure.h"
#include "tools/AdrOp.h"
#include "tools/constants.h"

#define FDTD_FLOAT float

class CartOperator
{
public:
    CartOperator();
    virtual ~CartOperator();

    void SetGeometryCSX(ContinuousStructure* geo);

	int CalcECOperator();

	void ApplyElectricBC(bool* dirs); //applied by default to all boundaries
	void ApplyMagneticBC(bool* dirs);

	double GetTimestep() {return dT;};

    void Reset();

protected:
    void Init();
    ContinuousStructure* CSX;

    AdrOp* MainOp;
    AdrOp* DualOp;
    double* discLines[3];
    unsigned int numLines[3];
	double gridDelta;

	double dT; //FDTD timestep!

	//EC operator
	FDTD_FLOAT* vv[3]; //calc new voltage from old voltage
	FDTD_FLOAT* vi[3]; //calc new voltage from old current
	FDTD_FLOAT* ii[3]; //calc new current from old current
	FDTD_FLOAT* iv[3]; //calc new current from old voltage

	//E-Field Excitation
	//!	  Calc the electric field excitation.
	bool CalcEFieldExcitation();
	unsigned int E_Ex_Count;
	unsigned int* E_Ex_index;
	FDTD_FLOAT* E_Ex_amp[3]; //represented as edge-voltages!!
	FDTD_FLOAT* E_Ex_delay;

	//Calc timestep only internal use
	double CalcTimestep();

	//EC elements, internal only!
    bool Calc_EC();
	bool Calc_ECPos(int n, unsigned int* pos, double* inEC);
	bool Calc_EffMatPos(int n, unsigned int* pos, double* inMat);
	double* EC_C[3];
    double* EC_G[3];
    double* EC_L[3];
    double* EC_R[3];
};

#endif // CARTOPERATOR_H
