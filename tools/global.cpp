/*
*	Copyright (C) 2024 Yifeng Li <tomli@tomli.me>
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
namespace po = boost::program_options;

// create global object
Global g_settings;

//! \brief This function initializes the object
Global::Global()
{
	m_showProbeDiscretization = false;
	m_nativeFieldDumps = false;
	m_VerboseLevel = 0;

}

po::options_description
Global::optionDesc()
{
	po::options_description optdesc("Additional global arguments");
	optdesc.add_options()
		(
			"showProbeDiscretization",
			po::bool_switch()->notifier(
				[&](bool val)
				{
					if (!val) return;
					cout << "openEMS - showing probe discretization information" << endl;
					m_showProbeDiscretization = true;
				}
			),
			"Show probe discretization information"
		)
		(
			"nativeFieldDumps",
			po::bool_switch()->notifier(
				[&](bool val)
				{
					if (!val) return;
					cout << "openEMS - dumping all fields using the native field components" << endl;
					m_nativeFieldDumps = true;
				}
			),
			"Dump all fields using the native field components"
		)
		(
			"verbose,v",
			po::value<unsigned int>()->default_value(0)->implicit_value(1)->
			notifier(
				[&](unsigned int val)
				{
					// Don't apply settings if the default value 0 is unchanged,
					// Apply settings and print messages if we have a non-default value
					// or if the non-default value is changed back to default in another
					// call (when running as a shared library).
					if (val == 0 && m_VerboseLevel == 0) return;

					m_VerboseLevel = val;
					cout << "openEMS - verbose level " << m_VerboseLevel << endl;
				}
			),
			"Verbose level, select debug level 1 to 3, "
			"also accept -v, -vv, -vvv"
		);
	return optdesc;
}

void Global::appendOptionDesc(po::options_description desc)
{
	m_optionDesc.add(desc);
}

void Global::parseLibraryArguments(std::vector<std::string> allOptions)
{
	clearOptions();

	for (std::string& option : allOptions)
	{
		if (option.length() == 1)
			option = "-" + option;
		else
			option = "--" + option;
	}

	// may throw
	po::store(
		po::command_line_parser(allOptions).options(m_optionDesc)
			.style(
				po::command_line_style::unix_style |
				po::command_line_style::case_insensitive)
			.run(),
		m_options
	);

	// run all registered callback functions in m_optionDesc
	po::notify(m_options);
}

void Global::parseCommandLineArguments(int argc, const char* argv[])
{
	// Hack: boost::program_options doesn't support repeated "-vv"
	// and "-vvv" syntax and causes validation failure. It's
	// not worthwhile to write a custom validator for exactly a
	// single special case. Just change argv[] to avoid them.
	std::pair<std::string, std::string> replaceTable[] =
	{
		{"-vv",  "--verbose=2"},
		{"-vvv", "--verbose=3"},
	};

	for (int i = 0; i < argc; i++)
	{
		for (const auto& entry : replaceTable)
		{
			if (std::string(argv[i]) == entry.first)
				argv[i] = entry.second.c_str();
		}
	}

	// may throw
	po::store(
		po::command_line_parser(argc, argv).options(m_optionDesc)
			.style(
				po::command_line_style::unix_style |
				po::command_line_style::case_insensitive)
			.run(),
		m_options
	);

	// run all registered callback functions in m_optionDesc
	po::notify(m_options);
}

void Global::showOptionUsage(std::ostream& ostr)
{
	ostr << m_optionDesc << endl;
}

bool Global::hasOption(std::string option)
{
	if (m_options.count(option) > 0)
		return true;
	else
		return false;
}

po::variable_value Global::getOption(std::string option)
{
	return m_options[option];
}

void Global::clearOptions()
{
	m_options.clear();
}
