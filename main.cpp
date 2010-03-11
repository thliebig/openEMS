#include <iostream>
#include <fstream>
#include <sstream>

#include "openems.h"

#include "examples/FDTD_examples.h"

#define STANDALONE

using namespace std;

int main(int argc, char *argv[])
{
	openEMS FDTD;

#ifdef STANDALONE
	if (argc<=1)
	{
		cerr << " usage: openEMS FDTD_XML_FILE [--disable-dumps]" << endl;
		exit(-1);
	}

	if (argc>=3)
	{
		for (int n=2;n<argc;++n)
		{
			if (strcmp(argv[n],"--disable-dumps")==0)
				FDTD.SetEnableDumps(false);
		}
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

	int EC = FDTD.SetupFDTD(file);
	if (EC) return EC;
	FDTD.RunFDTD();

	return 0;
}
