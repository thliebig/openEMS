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

class Global
{
public:
    Global();
	bool parseCommandLineArgument( const char *argv );

	bool showProbeDiscretization() const {return m_showProbeDiscretization;}

protected:
	bool m_showProbeDiscretization;
};

extern Global g_settings;

#endif // GLOBAL_H
