/*
 * Copyright (C) 2026 Thorsten Liebig <thorsten.liebig@gmx.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ARRAYLIB_MEMTRACK_H
#define ARRAYLIB_MEMTRACK_H

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <map>
#include <mutex>
#include <ostream>
#include <string>
#include <vector>

namespace ArrayLib
{
	class MemTrack;
};

// Accounting of the memory held by ArrayLib arrays: the total currently
// allocated, its high-water mark, and the per-array breakdown of what was
// live when that high-water mark was reached.
//
// This counts what the arrays asked for, not the resident set size. The
// allocator's own overhead, everything allocated outside ArrayLib (geometry,
// output buffers, thread stacks) and pages that were never touched are not
// included -- compare against PeakResidentBytes() for the whole picture.
class ArrayLib::MemTrack
{
public:
	struct Entry
	{
		size_t bytes;  // currently allocated by all arrays of this name
		size_t count;  // number of such arrays
		Entry() : bytes(0), count(0) {}
	};
	typedef std::map<std::string, Entry> Table;

	// Restores a stream's formatting on destruction, so that the report does
	// not leak fixed/precision into everything printed after it.
	class FormatGuard
	{
	public:
		FormatGuard(std::ostream& os) : m_os(os), m_flags(os.flags()), m_prec(os.precision()) {}
		~FormatGuard() {m_os.flags(m_flags); m_os.precision(m_prec);}
	protected:
		std::ostream& m_os;
		std::ios_base::fmtflags m_flags;
		std::streamsize m_prec;
	};

	static void Add(const std::string& name, size_t bytes)
	{
		State& s = GetState();
		std::lock_guard<std::mutex> lock(s.mtx);
		Entry& e = s.live[name];
		e.bytes += bytes;
		++e.count;
		s.current += bytes;
		if (s.current > s.peak)
		{
			s.peak = s.current;
			s.atPeak = s.live;
		}
	}

	static void Sub(const std::string& name, size_t bytes)
	{
		State& s = GetState();
		std::lock_guard<std::mutex> lock(s.mtx);
		Table::iterator it = s.live.find(name);
		if (it != s.live.end())
		{
			it->second.bytes -= bytes;
			if (--it->second.count == 0)
				s.live.erase(it);
		}
		s.current -= bytes;
	}

	//! Bytes currently held by all ArrayLib arrays
	static size_t Current()
	{
		State& s = GetState();
		std::lock_guard<std::mutex> lock(s.mtx);
		return s.current;
	}

	//! Highest value Current() has ever reached
	static size_t Peak()
	{
		State& s = GetState();
		std::lock_guard<std::mutex> lock(s.mtx);
		return s.peak;
	}

	//! Per-array breakdown of what is allocated now
	static Table Live()
	{
		State& s = GetState();
		std::lock_guard<std::mutex> lock(s.mtx);
		return s.live;
	}

	//! Per-array breakdown of what was allocated when Peak() was reached
	static Table AtPeak()
	{
		State& s = GetState();
		std::lock_guard<std::mutex> lock(s.mtx);
		return s.atPeak;
	}

	//! Print a table of arrays, largest first. numCells>0 adds a byte/cell column.
	static void PrintTable(std::ostream& os, const Table& table, size_t numCells=0)
	{
		FormatGuard guard(os);

		std::vector<Table::const_iterator> sorted;
		for (Table::const_iterator it=table.begin(); it!=table.end(); ++it)
			sorted.push_back(it);
		std::sort(sorted.begin(), sorted.end(), LargerFirst);

		os << "  " << std::left << std::setw(24) << "array" << std::right
		   << std::setw(6) << "count" << std::setw(14) << "MiB";
		if (numCells>0)
			os << std::setw(12) << "Byte/cell";
		os << std::endl;

		for (size_t n=0; n<sorted.size(); ++n)
		{
			const Entry& e = sorted.at(n)->second;
			os << "  " << std::left << std::setw(24) << sorted.at(n)->first << std::right
			   << std::setw(6) << e.count
			   << std::setw(14) << std::fixed << std::setprecision(2) << MiB(e.bytes);
			if (numCells>0)
				os << std::setw(12) << std::fixed << std::setprecision(2) << (double)e.bytes/(double)numCells;
			os << std::endl;
		}
	}

	//! Print the current and peak totals, with the breakdown of both
	static void Print(std::ostream& os, size_t numCells=0)
	{
		FormatGuard guard(os);

		Table live   = Live();
		Table atPeak = AtPeak();

		os << "Currently allocated\t: " << std::fixed << std::setprecision(2) << MiB(Current()) << " MiB" << std::endl;
		PrintTable(os, live, numCells);
		os << "Peak allocated\t\t: " << std::fixed << std::setprecision(2) << MiB(Peak()) << " MiB";
		if (numCells>0)
			os << " (" << std::fixed << std::setprecision(2) << (double)Peak()/(double)numCells << " Byte/cell)";
		os << std::endl;
		PrintTable(os, atPeak, numCells);
	}

protected:
	struct State
	{
		std::mutex mtx;
		size_t current;
		size_t peak;
		Table live;    // allocated right now
		Table atPeak;  // allocated when the high-water mark was reached
		State() : current(0), peak(0) {}
	};

	// Deliberately leaked: arrays may be freed by a destructor running after
	// a function-local static would already be gone.
	static State& GetState()
	{
		static State* s = new State();
		return *s;
	}

	static double MiB(size_t bytes) {return (double)bytes/1024.0/1024.0;}

	static bool LargerFirst(Table::const_iterator a, Table::const_iterator b)
	{
		return a->second.bytes > b->second.bytes;
	}
};

#endif // ARRAYLIB_MEMTRACK_H
