#include <iostream>
#include <time.h>
#include "FDTD/cartoperator.h"
#include "ContinuousStructure.h"

using namespace std;

int main(int argc, char *argv[])
{
    cout << "Hello World" << endl;
	fprintf(stderr,"%e\n",1.4456);

	time_t startTime=time(NULL);

    ContinuousStructure CSX;
	CSX.ReadFromXML("csx-files/1Mill.xml");

	CartOperator cop;
	cop.SetGeometryCSX(&CSX);

	cop.CalcECOperator();

	unsigned int* index = NULL;
	FDTD_FLOAT** amp=NULL;
	FDTD_FLOAT* delay=NULL;

	unsigned int nEx = cop.GetVoltageExcitation(index,amp,delay);

	time_t OpDoneTime=time(NULL);

	cerr << "Time for operator: " << difftime(OpDoneTime,startTime) << endl;

    return 0;
}
