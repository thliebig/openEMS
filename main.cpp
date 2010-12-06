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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include "openems.h"
#include "tools/global.h"

#ifndef GIT_VERSION
#define GIT_VERSION "unknown:compiled@" __DATE__
#endif

using namespace std;

int main(int argc, char *argv[])
{
	openEMS FDTD;

#ifdef _LP64
	string bits = "64bit";
#else
	string bits = "32bit";
#endif

	cout << " -------------------------------------------------------------------- " << endl;
	cout << " | openEMS " << bits << " -- version " GIT_VERSION << endl;
	cout << " | (C) 2010 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license" << endl;
	cout << " -------------------------------------------------------------------- " << endl;
	cout << openEMS::GetExtLibsInfo() << endl;

	if (argc<=1)
	{
		cout << " Usage: openEMS <FDTD_XML_FILE> [<options>...]" << endl << endl;
		cout << " <options>" << endl;
		cout << "\t--disable-dumps\t\tDisable all field dumps for faster simulation" << endl;
		cout << "\t--debug-material\tDump material distribution to a vtk file for debugging" << endl;
		cout << "\t--debug-operator\tDump operator to vtk file for debugging" << endl;
		cout << "\t--debug-boxes\t\tDump e.g. probe boxes to vtk file for debugging" << endl;
		cout << "\t--debug-CSX\t\tWrite CSX geometry file to debugCSX.xml" << endl;
		cout << "\t--engine=<type>\t\tChoose engine type" << endl;
		cout << "\t\t--engine=fastest\t\tfastest available engine" << endl;
		cout << "\t\t--engine=sse\t\t\tengine using sse vector extensions" << endl;
		cout << "\t\t--engine=sse_compressed\t\tengine using compressed operator + sse vector extensions" << endl;
		cout << "\t\t--engine=multithreaded\t\tengine using compressed operator + sse vector extensions + multithreading" << endl;
		cout << "\t--numThreads=<n>\tForce use n threads for multithreaded engine (needs: --engine=multithreaded)" << endl;
		cout << "\t--no-simulation\tonly run preprocessing; do not simulate" << endl;
		cout << "\n\t Additional global arguments " << endl;
		g_settings.ShowArguments(cout,"\t");
		cout << endl;
		exit(-1);
	}

	if (argc>=3)
	{
		for (int n=2; n<argc; ++n)
		{
			if (g_settings.parseCommandLineArgument(argv[n]))
				continue;

			if (!FDTD.parseCommandLineArgument(argv[n]))
				cout << "openEMS - unknown argument: " << argv[n] << endl;
		}
	}

	int EC = FDTD.SetupFDTD(argv[1]);
	if (EC) return EC;
	FDTD.RunFDTD();

	return 0;
}
