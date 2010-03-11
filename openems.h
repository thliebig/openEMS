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

protected:
	//! Number of Timesteps
	int NrTS;
	Operator* FDTD_Op;
	Engine* FDTD_Eng;
	ProcessingArray* PA;
};

#endif // OPENEMS_H
