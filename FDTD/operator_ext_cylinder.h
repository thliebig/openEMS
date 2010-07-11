#ifndef OPERATOR_EXT_CYLINDER_H
#define OPERATOR_EXT_CYLINDER_H

#include "operator_extension.h"
#include "operator.h"

class Operator_Cylinder;

class Operator_Ext_Cylinder : public Operator_Extension
{
	friend class Engine_Ext_Cylinder;
public:
	Operator_Ext_Cylinder(Operator_Cylinder* op);
	~Operator_Ext_Cylinder();

	virtual bool BuildExtension();

	virtual Engine_Extension* CreateEngineExtention();

	virtual bool IsCylinderCoordsSave() const {return true;}

	virtual std::string GetExtensionName() const {return std::string("Extension for the Cylinder-Coords Operator");}

	virtual void ShowStat(ostream &ostr) const;

protected:
	Operator_Cylinder* m_Op_Cyl;

	bool CC_closedAlpha;
	bool CC_R0_included;

	//special EC operator for R0
	FDTD_FLOAT* vv_R0; //calc new voltage from old voltage
	FDTD_FLOAT* vi_R0; //calc new voltage from old current

};

#endif // OPERATOR_EXT_CYLINDER_H
