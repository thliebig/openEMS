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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <sstream>
#include <boost/program_options.hpp>
#define _USE_MATH_DEFINES

#include "openems_global.h"

// declare a parameter as unused
#define UNUSED(x) (void)(x);

class OPENEMS_EXPORT Global
{
public:
	Global();

	bool showProbeDiscretization() const {return m_showProbeDiscretization;}

	//! Returns true if native field dumps are requested...
	bool NativeFieldDumps() const {return m_nativeFieldDumps;}
	//! Set dumps to use native fields.
	void SetNativeFieldDumps(bool val) {m_nativeFieldDumps=val;}

	//! Set the verbose level
	void SetVerboseLevel(int level) {m_VerboseLevel=level;m_SavedVerboseLevel=level;}
	//! Get the verbose level
	int GetVerboseLevel() const {return m_VerboseLevel;}

	//! Set a new verbose level temporarily, restore it with RestoreVerboseLevel()
	void SetTempVerboseLevel(int level) {m_SavedVerboseLevel=m_VerboseLevel;m_VerboseLevel=level;}
	//! Restore the temporarily overwritten verbose level
	void RestoreVerboseLevel() {m_VerboseLevel=m_SavedVerboseLevel;}

	// openEMS has different modules, and some options may be relevant
	// only to a single module, such as the upcoming Tiling engine.
	// Putting all getters and setters into the main openEMS class would
	// decreases modularity. It's the best if each module can handle their
	// own options. On the other hand, we really need to know all the options
	// globally to (1) check whether an option is valid, and (2) provide a
	// single API to set runtime options when running as a shared library
	// for Python binding.
	//
	// Thus, we use the following 3-step process.
	//
	// 1. If a class accepts options, it provides a method optionDesc()
	// that returns a options_description. During early initialization,
	// the main program openems.cpp collects them from multiple classes,
	// then registers them at here by calling appendOptionDesc(). Callback
	// functions can also be registered via options_description, which
	// may be used to set the internal state of a class.
	//
	// 2. Before starting simulation, if we're running as an executable,
	// parseCommandLineArguments() is called. If we're running as a shared
	// library, parseOption() is called (by the user). Callback functions
	// are also executed at this moment.
	//
	// 3. The parsed options is stored as a variables_map in m_options, which
	// can be accessed globally. Each module only needs to care about its
	// own options. If a class is only newed after parsing all options
	// (such as an engine-specific option), callbacks can't be used so the
	// getOption() method is used instead.

	// Return a list of supported options. If a class accepts options,
	// it should have its own method.
	boost::program_options::options_description optionDesc();

	// The optionDesc() of all classes are collected by openems.cpp by
	// calling appendOptionDesc.
	void appendOptionDesc(boost::program_options::options_description desc);

	// Parse all options provided as a std::vector of std::strings,
	// used when running as a shared library. May throw.
	void parseLibraryArguments(std::vector<std::string> allOptions);

	// Parse all options provided as C string in argv, used when running
	// as an executable. May throw.
	void parseCommandLineArguments(int argc, const char* argv[]);

	// Print usage of all known options
	void showOptionUsage(std::ostream& ostr);

	// Set, get, and clear options, provide access of options to all modules
	// globally.
	bool hasOption(std::string option);
	boost::program_options::variable_value getOption(std::string option);
	void clearOptions();

protected:
	bool m_showProbeDiscretization;
	bool m_nativeFieldDumps;
	int m_VerboseLevel;
	int m_SavedVerboseLevel;

	boost::program_options::variables_map m_options;
	boost::program_options::options_description m_optionDesc;
};

OPENEMS_EXPORT extern Global g_settings;

#endif // GLOBAL_H
