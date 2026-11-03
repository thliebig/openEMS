/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2026 openEMS contributors
 */
#include "openems.h"
#include "Common/processing.h"
#include "Common/processfields_fd.h"
#include "FDTD/engine.h"
#include "FDTD/operator.h"
#include "tinyxml.h"
#include <cstdio>
#include <fstream>
#include <set>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#endif

namespace
{
bool MakeDirectory(const std::string& path)
{
#ifdef _WIN32
	return _mkdir(path.c_str()) == 0;
#else
	return mkdir(path.c_str(), 0700) == 0;
#endif
}

std::string Precise(double value)
{
	std::ostringstream text;
	text << std::setprecision(17) << value;
	return text.str();
}

bool ValidID(const std::string& id)
{
	if (id.empty() || id.size() > 64)
	{
		return false;
	}
	for (char ch : id)
	{
		if (!((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
		      (ch >= '0' && ch <= '9') || ch == '-' || ch == '_'))
		{
			return false;
		}
	}
	return true;
}
}

void openEMS::CheckDFTSnapshot(double energy, double maxEnergy)
{
	std::ifstream request("DFT_SNAPSHOT", std::ios::binary);
	if (!request)
	{
		return;
	}
	std::string id;
	char ch;
	while (id.size() <= 65 && request.get(ch))
	{
		id += ch;
	}
	request.close();
	if (!id.empty() && id.back() == '\n')
	{
		id.pop_back();
	}
	if (!id.empty() && id.back() == '\r')
	{
		id.pop_back();
	}
	if (!ValidID(id))
	{
		std::cerr << "DFT snapshot: invalid request ID "
		          << "(use 1-64 ASCII letters/digits/-/_)." << std::endl;
		std::remove("DFT_SNAPSHOT");
		return;
	}
	// A single producer publishes requests atomically. Do not accept another while
	// this synchronous checkpoint is being written. Stale active files need review.
	std::ifstream active("DFT_SNAPSHOT.active");
	if (active || std::rename("DFT_SNAPSHOT", "DFT_SNAPSHOT.active") != 0)
	{
		return;
	}
	const auto started = std::chrono::steady_clock::now();
	const std::string root = "dft-snapshots";
	MakeDirectory(root);
	const std::string stage = root + "/." + id + ".tmp";
	const std::string target = root + "/" + id;
	struct stat info;
	bool ok = stat(target.c_str(), &info) != 0 && MakeDirectory(stage);
	TiXmlDocument manifest;
	TiXmlElement* header = new TiXmlElement("DFTSnapshot");
	header->SetAttribute("id", id);
	header->SetAttribute("build", GIT_VERSION);
	header->SetAttribute("timestep",
	                     std::to_string(FDTD_Eng->GetNumberOfTimesteps()));
	header->SetAttribute("time",
	                     Precise(FDTD_Eng->GetNumberOfTimesteps() * FDTD_Op->GetTimestep()));
	header->SetAttribute("lastEnergyEstimate", Precise(energy));
	header->SetAttribute("maxEnergyEstimate", Precise(maxEnergy));
	header->SetAttribute("endCriteria", Precise(endCrit));
	header->SetAttribute("converged", "false");
	manifest.LinkEndChild(header);
	std::set<std::string> filenames;
	if (ok)
	{
		// Drain every owner, including asynchronous TD dumps: HDF5 builds need not
		// be thread safe. No more samples are submitted until this returns.
		for (size_t n = 0; n < PA->GetNumberOfProcessings(); ++n)
		{
			PA->GetProcessing(n)->FinishAsync();
		}
		for (size_t n = 0; n < PA->GetNumberOfProcessings() && ok; ++n)
		{
			ProcessFieldsFD* fd = dynamic_cast<ProcessFieldsFD*>(PA->GetProcessing(n));
			if (!fd || !fd->SupportsDFTSnapshot())
			{
				continue;
			}
			std::string filename;
			ok = fd->WriteDFTSnapshot(stage, filename);
			ok = filenames.insert(filename).second && ok;
			TiXmlElement* file = new TiXmlElement("File");
			file->SetAttribute("name", filename);
			header->LinkEndChild(file);
		}
		ok = !filenames.empty() && ok;
		header->SetAttribute("copyWriteSeconds",
		                     Precise(std::chrono::duration<double>(
		                         std::chrono::steady_clock::now() - started).count()));
		ok = manifest.SaveFile((stage + "/snapshot.xml").c_str()) && ok;
		if (ok)
		{
			ok = std::rename(stage.c_str(), target.c_str()) == 0;
		}
	}
	if (ok)
	{
		std::cout << "DFT snapshot: published " << target << " at timestep "
		          << FDTD_Eng->GetNumberOfTimesteps() << "; continuing FDTD." << std::endl;
	}
	else
	{
		// Retain failed staging data for diagnosis; never publish a partial set.
		std::ofstream error(root + "/" + id + ".error");
		error << "DFT snapshot failed; inspect solver log and staging directory.\n";
		std::cerr << "DFT snapshot: failed " << id << "; continuing FDTD." << std::endl;
	}
	std::remove("DFT_SNAPSHOT.active");
}
