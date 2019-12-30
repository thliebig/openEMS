/*
*	Copyright (C) 2019 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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

#ifndef OPENEMS_GLOBAL_H
#define OPENEMS_GLOBAL_H

#if defined(WIN32)
	#ifdef BUILD_OPENEMS_LIB
	#define OPENEMS_EXPORT __declspec(dllexport)
	#else
	#define OPENEMS_EXPORT __declspec(dllimport)
	#endif
#else
#define OPENEMS_EXPORT
#endif

#endif // OPENEMS_GLOBAL_H
