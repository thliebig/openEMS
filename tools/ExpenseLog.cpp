#include "ExpenseLog.h"

ExpenseModule::ExpenseModule(const char* moduleName)
{
	chModuleName=moduleName;
	uiDoubleAdditions=uiDoubleMultiplications=uiIntAdditions=uiIntMultiplications=uiAssignments=uiBoolOp=0;
	uiMrdDA=uiMrdDM=uiMrdIA=uiMrdIM=uiMrdAssign=uiMrdBO=0;
}

ExpenseModule::~ExpenseModule() {};

void ExpenseModule::Clear()
{
	uiDoubleAdditions=uiDoubleMultiplications=uiIntAdditions=uiIntMultiplications=uiAssignments=uiBoolOp=0;
	uiMrdDA=uiMrdDM=uiMrdIA=uiMrdIM=uiMrdAssign=uiMrdBO=0;
}

void ExpenseModule::AddDoubleAdditons(unsigned int number)
{
	uiDoubleAdditions+=number;
	if (uiDoubleAdditions>=MRD)
	{
		uiDoubleAdditions-=MRD;
		++uiMrdDA;
	}
}

void ExpenseModule::AddDoubleMultiplications(unsigned int number)
{
	uiDoubleMultiplications+=number;
	if (uiDoubleMultiplications>=MRD)
	{
		uiDoubleMultiplications-=MRD;
		++uiMrdDM;
	}
}

void ExpenseModule::AddIntAdditons(unsigned int number)
{
	uiIntAdditions+=number;
	if (uiIntAdditions>=MRD)
	{
		uiIntAdditions-=MRD;
		++uiMrdIA;
	}
}

void ExpenseModule::AddIntMultiplications(unsigned int number)
{
	uiIntMultiplications+=number;
	if (uiIntMultiplications>=MRD)
	{
		uiIntMultiplications-=MRD;
		++uiMrdIM;
	}
}

void ExpenseModule::AddAssignments(unsigned int number)
{
	uiAssignments+=number;
	if (uiAssignments>=MRD)
	{
		uiAssignments-=MRD;
		++uiMrdAssign;
	}
}

void ExpenseModule::AddBoolOperations(unsigned int number)
{
	uiBoolOp+=number;
	if (uiBoolOp>=MRD)
	{
		uiBoolOp-=MRD;
		++uiMrdBO;
	}
}

void ExpenseModule::AddOperations(unsigned int IntAdd, unsigned int IntMul, unsigned int DoubleAdd, unsigned int DoubleMul, unsigned int Assigns, unsigned int BoolOp)
{
	this->AddIntAdditons(IntAdd);
	this->AddIntMultiplications(IntMul);
	this->AddDoubleAdditons(DoubleAdd);
	this->AddDoubleMultiplications(DoubleMul);
	this->AddAssignments(Assigns);
	this->AddBoolOperations(BoolOp);
}

void ExpenseModule::PrintfSelf(FILE* file)
{
	fprintf(file,"\n***********\n Module:  %s\n Additions:\n Double: %3.0d%9d\tInteger: %3.0d%9d",chModuleName,uiMrdDA,uiDoubleAdditions,uiMrdIA,uiIntAdditions);
	fprintf(file,"\n\n Multiplications:\n Double: %3.0d%9d\tInteger: %3.0d%9d\n",uiMrdDM,uiDoubleMultiplications,uiMrdIM,uiIntMultiplications);
	fprintf(file,"\n Assignments: %3.0d%9d\tBool Operations: %3.0d%9d\n",uiMrdAssign,uiAssignments,uiMrdBO,uiBoolOp);
	fprintf(file,"\n***********\n");
}




/***********************************************************************************************************************/

ExpenseLog::ExpenseLog(void)
{
}

ExpenseLog::~ExpenseLog(void)
{
	for (size_t i=0; i<vModules.size(); ++i)
	{
		delete vModules.at(i);
		vModules.at(i)=NULL;
	}
	vModules.clear();
}

ExpenseModule* ExpenseLog::AddModule(const char* name)
{
	ExpenseModule* newModule = new ExpenseModule(name);
	vModules.push_back(newModule);
	return newModule;
}

void ExpenseLog::PrintAll(FILE *file)
{
	double totalAdd=0,totalMul=0,totalBool=0,totalAssign=0;
	fprintf(stderr,"\n ----------------\n Expense Log PrintOut\n Nr of Modules: %d\n",vModules.size());
	for (size_t i=0; i<vModules.size(); ++i)
	{
		totalAdd+=((double)vModules.at(i)->uiIntAdditions+(double)vModules.at(i)->uiDoubleAdditions) + 1e9*((double)vModules.at(i)->uiMrdIA+(double)vModules.at(i)->uiMrdIA);
		totalMul+=((double)vModules.at(i)->uiIntMultiplications+(double)vModules.at(i)->uiDoubleMultiplications) + 1e9*(double)(vModules.at(i)->uiMrdIM+vModules.at(i)->uiMrdIM);
		totalBool+=(double)vModules.at(i)->uiBoolOp + 1e9*(double)vModules.at(i)->uiMrdBO;
		totalAssign+=(double)vModules.at(i)->uiAssignments + 1e9*(double)vModules.at(i)->uiMrdAssign;
		vModules.at(i)->PrintfSelf(file);
	}
	fprintf(stderr," Total:\n Additions: %e Multiplications: %e\n Bool Operations: %e Assignments: %e\n",totalAdd,totalMul,totalBool,totalAssign);
	fprintf(stderr,"\n ----------------\n");
}

void ExpenseLog::ClearAll()
{
	for (size_t i=0; i<vModules.size(); ++i)
	{
		vModules.at(i)->Clear();
	}
}

