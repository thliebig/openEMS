/*
*	Copyright (C) 2012-2014 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#include "nf2ff.h"

using namespace std;

int main(int argc, char *argv[])
{
	cout << " ---------------------------------------------------------------------- " << endl;
	cout << " | nf2ff, near-field to far-field transformation for openEMS " << endl;
	cout << " | (C) 2012-2014 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license" << endl;
	cout << " ---------------------------------------------------------------------- " << endl;

	if (argc<=1)
	{
		cout << " Usage: nf2ff <nf2ff-xml-file>" << endl << endl;
		cout << endl;
		exit(-1);
	}

	if (argc>=2)
	{
		return !nf2ff::AnalyseXMLFile(argv[argc-1]);
	}

	return 0;
}
