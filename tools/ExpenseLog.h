#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

#define EXPENSE_LOG 0
#define MRD 1000000000

#if EXPENSE_LOG==1

#define EXPENSE_DEFINE \
ExpenseLog EL; \
ExpenseModule* EngineExpense=EL.AddModule("Static Engine Expenses"); \
ExpenseModule* PPExpense=EL.AddModule("Static Post Processing"); \
ExpenseModule* AdrOpExpense=EL.AddModule("Adress Operator");
#define EXTERN_EXPENSE_DEFINE extern ExpenseLog EL;
#define ENGINEEXPENSE_DEFINE extern ExpenseModule* EngineExpense;
#define POSTPROCEXPENSE_DEFINE extern ExpenseModule* PPExpense;
#define ADREXPENSE_DEFINE extern ExpenseModule* AdrOpExpense;
#define ENGINEEXPENSE(IA,IM,DA,DM,AS,BO) EngineExpense->AddOperations((IA),(IM),(DA),(DM),(AS),(BO));
#define POSTPROCEXPENSE(IA,IM,DA,DM,AS,BO) PPExpense->AddOperations((IA),(IM),(DA),(DM),(AS),(BO));
#define ADRESSEXPENSE(IA,IM,DA,DM,AS,BO) AdrOpExpense->AddOperations((IA),(IM),(DA),(DM),(AS),(BO));
#define EXPENSEPRINT EL.PrintAll(stderr);
#define EXPENSECLEAR EL.ClearAll();
#else

#define EXPENSE_DEFINE
#define EXTERN_EXPENSE_DEFINE
#define ENGINEEXPENSE_DEFINE
#define POSTPROCEXPENSE_DEFINE
#define ADREXPENSE_DEFINE
#define ENGINEEXPENSE(IA,IM,DA,DM,AS,BO)
#define POSTPROCEXPENSE(IA,IM,DA,DM,AS,BO)
#define ADRESSEXPENSE(IA,IM,DA,DM,AS,BO)
#define EXPENSEPRINT
#define EXPENSECLEAR
#endif

class ExpenseModule
{
	friend class ExpenseLog;
public:
	ExpenseModule(const char* moduleName);
	~ExpenseModule();

	void Clear();

	void AddDoubleAdditons(unsigned int number);
	void AddDoubleMultiplications(unsigned int number);

	void AddIntAdditons(unsigned int number);
	void AddIntMultiplications(unsigned int number);

	void AddAssignments(unsigned int number);
	void AddBoolOperations(unsigned int number);

	void AddOperations(unsigned int IntAdd, unsigned int IntMul, unsigned int DoubleAdd, unsigned int DoubleMul, unsigned int Assigns, unsigned int BoolOp);

	void PrintfSelf(FILE* file=stdout);

protected:
	const char* chModuleName;
	unsigned int uiDoubleAdditions;
	unsigned int uiDoubleMultiplications;
	unsigned int uiIntAdditions;
	unsigned int uiIntMultiplications;
	unsigned int uiAssignments;
	unsigned int uiBoolOp;
	unsigned int uiMrdDA;
	unsigned int uiMrdDM;
	unsigned int uiMrdIA;
	unsigned int uiMrdIM;
	unsigned int uiMrdAssign;
	unsigned int uiMrdBO;
};

class ExpenseLog
{
public:
	ExpenseLog(void);
	~ExpenseLog(void);

	ExpenseModule* AddModule(const char* name);
	void PrintAll(FILE *file=stdout);
	void ClearAll();
protected:
	vector<ExpenseModule*> vModules;
};


