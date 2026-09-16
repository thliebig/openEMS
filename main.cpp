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

#include "tools/global.h"

#ifndef GIT_VERSION
#define GIT_VERSION "unknown:compiled@" __DATE__
#endif

using namespace std;

int main(int argc, const char* argv[])
{
	openEMS FDTD;

	openEMS::WelcomeScreen();

	if (argc<=1)
	{
		FDTD.showUsage();
		return -1;
	}

	g_settings.parseCommandLineArguments(argc, argv);

	int EC = FDTD.ParseFDTDSetup(argv[1]);
	if (!EC) {
		cerr << "openEMS - ParseFDTDSetup failed." << endl;
		return 1;
	}
	EC = FDTD.SetupFDTD();
	if (EC) return EC;
	FDTD.RunFDTD();

	return 0;
}
