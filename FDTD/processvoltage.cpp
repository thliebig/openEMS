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
	FDTD_FLOAT voltage=0;
	unsigned int pos[3]={start[0],start[1],start[2]};
//	cerr << Eng->volt[1][pos[0]][pos[1]][pos[2]] << endl;
	for (int n=0;n<3;++n)
	{
		if (start[n]<stop[n])
		{
			for (;pos[n]<stop[n];++pos[n])
			{
				voltage+=Eng->volt[n][pos[0]][pos[1]][pos[2]];
//				cerr << n << " " << pos[0] << " " << pos[1] << " " <<  pos[2] << " " << Eng->volt[n][pos[0]][pos[1]][pos[2]] << endl;
			}

		}
		else
		{
			for (;pos[n]>stop[n];--pos[n])
				voltage-=Eng->volt[n][pos[0]][pos[1]][pos[2]];
		}
	}
//	cerr << voltage << endl;
	voltages.push_back(voltage);
	file << (double)Eng->numTS*Op->GetTimestep() << "\t" << voltage << endl;
	return GetNextInterval();
}
