/*
*	Copyright (C) 2025 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include <boost/program_options.hpp>
using namespace boost::program_options;

#include "sar_calculation.h"

using namespace std;

int main(int argc, char *argv[])
{
	cout << " ---------------------------------------------------------------------- " << endl;
	cout << " | SAR calculation for openEMS "                                          << endl;
	cout << " | (C) 2012-2025 Thorsten Liebig <thorsten.liebig@gmx.de>  GPL license"   << endl;
	cout << " ---------------------------------------------------------------------- " << endl;

	options_description desc("Options");
	variables_map vm;
	try {
		string ifile, ofile, method;
		double m_avg = 0;
		bool debug = false;
		bool export_cube_stats = false;
		bool legacyHDF5 = false;

		desc.add_options()
		("help,h"  , "print usage message")
		("input,i" , value(&ifile)->required(), "pathname to input hdf5 file")
		("output,o", value(&ofile)->required(), "pathname for output hdf5 file")
		("method",   value(&method)->default_value("SIMPLE"), "set SAR method: IEEE_C95_3, IEEE_62704, SIMPLE")
		("mass,m"  , value(&m_avg), "averaging mass in g")
		("verbose,v", bool_switch(&debug), "verbose")
		("export_cube_stats,e", bool_switch(&export_cube_stats), "Export Cube Statistics")
		("legacyHDF5Dumps", bool_switch(&legacyHDF5), "Dumping using the legacy HDF5 file format as required for Octave/Matlab import");

		store(parse_command_line(argc, argv, desc), vm);

		if (vm.count("help")) {
			cout << desc << "\n";
			return 0;
		}

		notify(vm);

		SAR_Calculation sar_calc;
		sar_calc.SetDebugLevel(int(debug));
		sar_calc.SetAveragingMass(m_avg/1000);
		//sar_calc.SetAveragingMethod(SAR_Calculation::IEEE_62704);
		if (!sar_calc.SetAveragingMethod(method, !debug))
			return -1;
		sar_calc.CalcFromHDF5(ifile, ofile, export_cube_stats, legacyHDF5);
		return 0;
	}
	catch(exception& e) {
		cerr << e.what() << "\n";
	}
}
