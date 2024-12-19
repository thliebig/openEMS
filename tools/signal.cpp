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

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include "signal.h"

static volatile std::sig_atomic_t m_sigintAbort = 0;

#ifndef WIN32
static void (*m_sigHandlerOriginal)(int) = NULL;
#else
static PHANDLER_ROUTINE m_sigHandlerRegistered = NULL;
#endif

void Signal::SetupHandlerForSIGINT(int type)
{
	m_sigintAbort = 0;

#ifndef WIN32
	UnixSetupHandlerForSIGINT(type);
#else
	Win32SetupHandlerForConsoleCtrl(type);
#endif
}

#ifndef WIN32
void Signal::UnixSetupHandlerForSIGINT(int type)
{
	if (type == SIGNAL_ORIGINAL && m_sigHandlerOriginal)
	{
		// If we're acting as a shared library and a program (such as
		// the Python interpreter) calls us, changing the SIGINT handler
		// unilaterally may overwrite the original handler and affect
		// the functionality of the original program. Thus, we save the
		// original handler and restore it after the end of SetupFDTD()
		// or RunFDTD() to minimize the disruption.
		auto retval = std::signal(SIGINT, m_sigHandlerOriginal);
		if (retval == SIG_ERR)
		{
			fprintf(stderr, "Signal::UnixSetupHandlerForSIGINT(): "
					"Failed to restore signal handler!\n");
		}
		m_sigHandlerOriginal = NULL;
	}
	else if (type == SIGNAL_EXIT_GRACEFUL)
	{
		m_sigHandlerOriginal = std::signal(SIGINT, UnixGracefulExitHandler);
		if (m_sigHandlerOriginal == SIG_ERR)
		{
			fprintf(stderr, "Signal::UnixSetupHandlerForSIGINT(): "
					"Failed to set UnixGracefulExitHandler!\n");
			m_sigHandlerOriginal = NULL;
		}
	}
	else if (type == SIGNAL_EXIT_FORCE)
	{
		m_sigHandlerOriginal = std::signal(SIGINT, UnixForceExitHandler);
		if (m_sigHandlerOriginal == SIG_ERR)
		{
			fprintf(stderr, "Signal::UnixSetupHandlerForSIGINT(): "
					"Failed to set UnixForceExitHandler!\n");
			m_sigHandlerOriginal = NULL;
		}
	}
}

void Signal::UnixGracefulExitHandler(int signal)
{
	m_sigintAbort = 1;

	// C standard only guarantees that a sig_atomic_t variable is safe
	// to read or write, but it's not necessarily safe to increment by
	// one, and also not safe to set one sig_atomic_t depending on the
	// result of another sig_atomic_t.
	//
	// Thus, we switch the signal handler itself instead of recording
	// the number of times SIGINT is raised.
	auto retval = std::signal(SIGINT, UnixForceExitHandler);
	if (retval == SIG_ERR)
	{
		SafeStderrWrite("\nSignal::UnixGracefulExitHandler(): "
				"Failed to set UnixForceExitHandler!");
	}
	else
	{
		SafeStderrWrite("\nSignal::UnixGracefulExitHandler(): "
				"Gracefully aborting simulation "
				"now, this may take a few seconds...\n"
				"Signal::UnixGracefulExitHandler(): "
				"To force-exit, send Ctrl-C again, "
				"but simulation results may be lost.\n");
	}
}

void Signal::UnixForceExitHandler(int signal)
{
	SafeStderrWrite("\nSignal::UnixForceExitHandler(): "
			"Force-exit simulation process now!\n");

	// By convention, if a program is (uncleanly) aborted due to
	// an external signal, preferably it should return 128 + signal.
	// For SIGINT, it's 130.
	std::_Exit(128 + signal);
}

#else

void Signal::Win32SetupHandlerForConsoleCtrl(int type)
{
	if (type == SIGNAL_ORIGINAL || m_sigHandlerRegistered)
	{
		// On Windows, SetConsoleCtrlHandler appends a new ConsoleCtrlHandler
		// in addition to the existing handlers. Thus, we need to record
		// the ConsoleCtrlHandler installed by us (instead of getting the
		// pre-existing handlers on Unix). Then, before we install a new
		// signal handler, we need to use the argument "Add == FALSE" to
		// remove the handler we previously installed.
		//
		// We also need to do the same in case that we're restoring the
		// ConsoleCtrlHandler to the original state (note how on Unix, the
		// if expression uses "AND", but on Windows, the if expression uses
		// "OR".
		BOOL success = SetConsoleCtrlHandler(m_sigHandlerRegistered, FALSE);
		m_sigHandlerRegistered = NULL;

		if (!success)
		{
			fprintf(stderr, "Signal::Win32SetupHandlerForConsoleCtrl(): "
					"Failed to unregister ConsoleCtrlHandler!\n");
			return;
		}
	}

	// Assume m_sigHandlerRegistered has already been unregistered.
	if (type == SIGNAL_EXIT_GRACEFUL)
	{
		m_sigHandlerRegistered = (PHANDLER_ROUTINE) Win32GracefulExitHandler;
		BOOL success = SetConsoleCtrlHandler(m_sigHandlerRegistered, TRUE);

		if (!success)
		{
			fprintf(stderr, "Signal::Win32SetupHandlerForConsoleCtrl(): "
					"Failed to register Win32GracefulExitHandler!\n");
		}
	}
	else if (type == SIGNAL_EXIT_FORCE)
	{
		m_sigHandlerRegistered = (PHANDLER_ROUTINE) Win32ForceExitHandler;
		BOOL success = SetConsoleCtrlHandler(m_sigHandlerRegistered, TRUE);

		if (!success)
		{
			fprintf(stderr, "Signal::Win32SetupHandlerForConsoleCtrl(): "
					"Failed to register Win32ForceExitHandler!\n");
		}
	}
}

BOOL Signal::Win32GracefulExitHandler(DWORD fdwCtrlType)
{
	m_sigintAbort = 1;

	// unregister the current handler
	BOOL success = SetConsoleCtrlHandler(m_sigHandlerRegistered, FALSE);
	if (!success)
	{
		SafeStderrWrite("Signal::Win32GracefulExitHandler(): "
				"Failed to unregister Win32GracefulExitHandler!\n");
		return true;
	}

	// install a new handler
	m_sigHandlerRegistered = (PHANDLER_ROUTINE) Win32ForceExitHandler;
	success = SetConsoleCtrlHandler(m_sigHandlerRegistered, TRUE);
	if (!success)
	{
		SafeStderrWrite("Signal::Win32GracefulExitHandler(): "
				"Failed to register Win32ForceExitHandler!\n");
	}
	else
	{
		SafeStderrWrite("\nSignal::Win32GracefulExitHandler(): "
				"Gracefully aborting simulation "
				"now, this may take a few seconds...\n"
				"Signal::Win32GracefulExitHandler(): "
				"To force-exit, send Ctrl-C again, "
				"but simulation results may be lost.\n");
	}

	return true;
}

BOOL Signal::Win32ForceExitHandler(DWORD fdwCtrlType)
{
	SafeStderrWrite("\nSignal::Win32ForceExitHandler(): "
			"Force-exit simulation process now!\n");

	// On Windows, the exit code for SIGINT is always 3.
	std::_Exit(3);

	// unreachable
	return true;
}
#endif

bool Signal::ReceivedSIGINT(void)
{
	if (m_sigintAbort)
		return true;
	else
		return false;
}

void Signal::SafeStderrWrite(const char *buf)
{
#ifdef WIN32
	// On Windows, using any kind of system calls in a ANSI C signal
	// handler is prohibited, in this case, this function should return
	// immediately without doing anything. But, when the official way
	// SetConsoleCtrlHandler() is used (instead of using ANSI C signals),
	// there's no such restriction.
	fprintf(stderr, "%s", buf);
	fflush(stderr);
	return;
#else
	// On Unix, in a signal handler, it's unsafe to use normal I/O
	// functions such as iostream, puts(), printf(), fprintf(). The
	// only safe option is the system call write().
	size_t buf_len = strlen(buf);
	ssize_t bytes = 0;

	while (buf_len > 0)
	{
		bytes = write(STDERR_FILENO, buf, buf_len);
		if (bytes < 0)
		{
			// write failure, nothing we can do.
			return;
		}

		if ((size_t) bytes > buf_len)
		{
			// Assertion: This should never happen. bytes is
			// always less or equal to buf_len, and buf_len
			// will never underflow under any circumstances
			// (unless the write system call is broken).
			return;
		}

		buf += bytes;		    // advance buffer position
		buf_len -= (size_t) bytes;  // decrement limiter
	}
	return;  // write completed.
#endif
}
