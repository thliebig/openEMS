/*
*	Copyright (C) 2010 Thorsten Liebig (Thorsten.Liebig@gmx.de)
*
*	This program is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*
*	This program is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*
*	You should have received a copy of the GNU General Public License
*	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <sstream>

#include "openems.h"

#include "examples/FDTD_examples.h"

#define STANDALONE

#ifndef GIT_VERSION
#define GIT_VERSION "unknown:compiled@" __DATE__
#endif

using namespace std;

int main(int argc, char *argv[])
{
	openEMS FDTD;

#ifdef STANDALONE
	if (argc<=1)
	{
		cout << " | openEMS version " GIT_VERSION << endl;
		cout << " | (C) 2010 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license" << endl << endl;
		cout << " usage: openEMS FDTD_XML_FILE [--disable-dumps] [--debug-material] [--debug-operator]" << endl;
		exit(-1);
	}

	if (argc>=3)
	{
		for (int n=2;n<argc;++n)
		{
			if (!FDTD.parseCommandLineArgument(argv[n]))
				cout << "openEMS - unknown argument: " << argv[n] << endl;
		}
	}

	char* file = argv[1];
#else
	//*************** setup/read geometry ************//

	const char* fileDP="../examples/Dipol.xml";
	BuildDipol(fileDP);

	const char* filePW="../examples/PlaneWave.xml";
	BuildPlaneWave(filePW);

	const char* fileMSL="../examples/MSL.xml";
	BuildMSL(fileMSL);

	const char* fileCoax="../examples/Coax_Cart.xml";
	BuildCoaxial_Cartesian(fileCoax);

	const char* fileHelix="../examples/Helix.xml";
	BuildHelix(fileHelix);

	const char* file=fileHelix;

	FDTD.DebugMaterial();
//	FDTD.DebugOperator();
#endif

	int EC = FDTD.SetupFDTD(file);
	if (EC) return EC;
	FDTD.RunFDTD();

	return 0;
}
