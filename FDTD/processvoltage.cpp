#include "processvoltage.h"

ProcessVoltage::ProcessVoltage(Operator* op, Engine* eng) : Processing(op, eng)
{
}

ProcessVoltage::~ProcessVoltage()
{
	file.close();
}

void ProcessVoltage::OpenFile(string outfile)
{
	if (file.is_open()) file.close();

	file.open(outfile.c_str());
	if (file.is_open()==false)
	{
		cerr << "Can't open file: " << outfile << endl;
		return;
	}
}

int ProcessVoltage::Process()
{
	if (Enabled==false) return -1;
	if (CheckTimestep()==false) return GetNextInterval();
	FDTD_FLOAT voltage=CalcLineIntegral(start,stop,0);
//	cerr << voltage << endl;
	voltages.push_back(voltage);
	file << (double)Eng->numTS*Op->GetTimestep() << "\t" << voltage << endl;
	return GetNextInterval();
}
