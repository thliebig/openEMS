#include "processing.h"

Processing::Processing(Operator* op, Engine* eng)
{
	Op=op;
	Eng=eng;
}

Processing::~Processing()
{
	file.close();
}

void Processing::OpenFile(string outfile)
{
	if (file.is_open()) file.close();

	file.open(outfile.c_str());
	if (file.is_open()==false)
	{
		cerr << "Can't open file: " << outfile << endl;
		return;
	}
}
