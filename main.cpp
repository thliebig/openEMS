#include <iostream>
#include <fstream>
#include <sstream>

#include "openems.h"

#include "examples/FDTD_examples.h"

//#define STANDALONE

using namespace std;

int main(int argc, char *argv[])
{
#ifdef STANDALONE
	if (argc<=1)
	{
		cerr << " argc= " << argc << endl;
		exit(-1);
	}

	char* file = argv[1];
#else
	//*************** setup/read geometry ************//

	const char* fileDP="examples/Dipol.xml";
	BuildDipol(fileDP);

	const char* filePW="examples/PlaneWave.xml";
	BuildPlaneWave(filePW);

	const char* fileMSL="examples/MSL.xml";
	BuildMSL(fileMSL);

	const char* fileCoax="examples/Coax_Cart.xml";
	BuildCoaxial_Cartesian(fileCoax);

	const char* file=fileCoax;

//	cerr << CSX.ReadFromXML("examples/PlaneWave.xml") << endl;
#endif
	openEMS FDTD;

	int EC = FDTD.SetupFDTD(file);
	if (EC) return EC;
	FDTD.RunFDTD();

	return 0;
}
