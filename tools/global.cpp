/*
*	Copyright (C) 2010 Sebastian Held <sebastian.held@gmx.de>
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

#include <cstring>
#include <iostream>
#include "global.h"

using namespace std;

// create global object
Global g_settings;

Global::Global()
{
	m_showProbeDiscretization = false;
}

//! \brief This function initializes the object
bool Global::parseCommandLineArgument( const char *argv )
{
	if (!argv)
		return false;

	if (strcmp(argv,"--showProbeDiscretization")==0)
	{
		cout << "openEMS - showing probe discretization information" << endl;
		m_showProbeDiscretization = true;
		return true;
	}

	return false;
}
