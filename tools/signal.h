/*
*	Copyright (C) 2023 Yifeng Li <tomli@tomli.me>
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

#ifndef SIGNAL_H
#define SIGNAL_H

#include <csignal>

#ifndef WIN32
#include <unistd.h>
#else
#include <windows.h>
#endif

enum
{
	SIGNAL_ORIGINAL,
	SIGNAL_EXIT_GRACEFUL,
	SIGNAL_EXIT_FORCE,
};

class Signal
{
public:
	static void SetupHandlerForSIGINT(int type);
	static bool ReceivedSIGINT(void);

private:
	static void SafeStderrWrite(const char *buf);

#ifndef WIN32
	static void UnixSetupHandlerForSIGINT(int type);
	static void UnixGracefulExitHandler(int signal);
	static void UnixForceExitHandler(int signal);
#else
	static void Win32SetupHandlerForConsoleCtrl(int type);
	static BOOL Win32GracefulExitHandler(DWORD fdwCtrlType);
	static BOOL Win32ForceExitHandler(DWORD fdwCtrlType);
#endif
};

#endif
