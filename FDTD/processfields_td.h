#ifndef PROCESSFIELDS_TD_H
#define PROCESSFIELDS_TD_H

#include "processfields.h"

class ProcessFieldsTD : public ProcessFields
{
public:
	ProcessFieldsTD(Operator* op, Engine* eng);
	virtual ~ProcessFieldsTD();

	virtual int Process();

	//! Set the length of the filename timestep pad filled with zeros (default is 8)
	void SetPadLength(int val) {pad_length=val;};

protected:
	int pad_length;

	void DumpNoInterpol(ofstream &file);
	void DumpCellInterpol(ofstream &file);
};

#endif // PROCESSFIELDS_TD_H
