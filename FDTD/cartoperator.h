#ifndef CARTOPERATOR_H
#define CARTOPERATOR_H

#include "operator.h"

class CartOperator : public Operator
{
public:
    CartOperator();
    virtual ~CartOperator();

	virtual void SetGeometryCSX(ContinuousStructure* geo);

	virtual int CalcECOperator();

	virtual void ApplyElectricBC(bool* dirs); //applied by default to all boundaries
	virtual void ApplyMagneticBC(bool* dirs);

	virtual void Reset();

protected:
	virtual void Init();

    AdrOp* MainOp;
    AdrOp* DualOp;

	virtual bool CalcEFieldExcitation();
	virtual bool CalcPEC();
	virtual double CalcTimestep();

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
