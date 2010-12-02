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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <sstream>

#define UNUSED(x) (void)(x);

class Global
{
public:
    Global();
	bool parseCommandLineArgument( const char *argv );

	bool showProbeDiscretization() const {return m_showProbeDiscretization;}

	//! Returns true if native field dumps are requested...
	bool NativeFieldDumps() const {return m_nativeFieldDumps;}
	//! Set dumps to use native fields.
	void SetNativeFieldDumps(bool val) {m_nativeFieldDumps=val;}

	void ShowArguments(std::ostream& ostr, std::string front=std::string());

	void SetVerboseLevel(int level) {m_VerboseLevel=level;}
	int GetVerboseLevel() const {return m_VerboseLevel;}

protected:
	bool m_showProbeDiscretization;
	bool m_nativeFieldDumps;
	int m_VerboseLevel;
};

extern Global g_settings;

// declare a parameter as unused
#define UNUSED(x) (void)(x);

#endif // GLOBAL_H
