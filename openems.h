#ifndef OPENEMS_H
#define OPENEMS_H

class Operator;
class Engine;
class ProcessingArray;

class openEMS
{
public:
    openEMS();
	~openEMS();

	int SetupFDTD(const char* file);

	void RunFDTD();

	void Reset();

	void SetEnableDumps(bool val) {Enable_Dumps=val;}

protected:
	//! Number of Timesteps
	int NrTS;
	bool Enable_Dumps;
	Operator* FDTD_Op;
	Engine* FDTD_Eng;
	ProcessingArray* PA;
};

#endif // OPENEMS_H
