/*
* Copyright (C) 2026 openEMS contributors
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#include "operator_metal.h"
#include "engine_metal.h"
#include "extensions/operator_ext_upml.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/operator_ext_mur_abc.h"
#include "tools/useful.h"
#include <algorithm>
#include <cstdint>
#include <limits>
#include <string>
#include <thread>
#include <typeinfo>

using std::cout;
using std::endl;

Operator_Metal* Operator_Metal::New(unsigned int threads)
{
	// Fail before the costly material/geometry sampling if no GPU is present.
	std::string reason;
	if (!MetalDeviceAvailable(reason))
	{
		std::cerr << "Metal: " << reason << "; cannot run the Metal engine" << std::endl;
		return nullptr;
	}
	cout << "Create FDTD operator (Metal field updates)" << endl;
	Operator_Metal* op = new Operator_Metal();
	op->Init();
	op->m_setupThreads = threads;
	return op;
}

Operator_Metal::Operator_Metal() : Operator_sse(), m_setupThreads(0),
	m_geoWinnersValid(false), m_geoWinnersWarned(false)
{
}

unsigned int Operator_Metal::GetSetupThreads() const
{
	return m_setupThreads ? m_setupThreads : std::max(1U,std::thread::hardware_concurrency());
}

const std::vector<Operator::GeometryWinner>* Operator_Metal::GetGeometryWinners(GeometryWinnerType type, bool dualMesh) const
{
	if (!m_geoWinnersValid)
	{
		if (!m_geoWinnersWarned)
		{
			std::cerr << "Metal: geometry winners unavailable (PEC pass was skipped or fell back to the CPU); "
			             "EC-consuming extensions fall back to full-grid CSXCAD lookups" << std::endl;
			m_geoWinnersWarned = true;
		}
		return nullptr;
	}
	switch (type)
	{
	case GEO_CONDUCTING_SHEET:
		return dualMesh ? nullptr : &m_geoConductingSheet;
	case GEO_DISPERSIVE:
		return dualMesh ? &m_geoDispersiveDual : &m_geoDispersivePrimal;
	}
	return nullptr;
}

bool Operator_Metal::SetupCSXGrid(CSRectGrid* grid)
{
	if (!Operator_sse::SetupCSXGrid(grid))
		return false;
	// The kernels address the field as packed float4 words (diamond update) and as
	// scalar components (excitation, ADE, indexed UPML). Reject a grid the 32-bit
	// index format cannot cover before the material and coefficient build.
	const uint64_t zSlots  = ((uint64_t)numLines[2] + 3) / 4;
	const uint64_t scalars = 12 * (uint64_t)numLines[0] * numLines[1] * zSlots;
	if (scalars > std::numeric_limits<uint32_t>::max())
	{
		std::cerr << "Metal: grid has "
		          << (uint64_t)numLines[0] * numLines[1] * numLines[2]
		          << " cells; the packed field index (" << scalars
		          << ") exceeds the 32-bit kernel format. Aborting before setup."
		          << std::endl;
		return false;
	}
	return true;
}

bool Operator_Metal::Calc_EC()
{
	if (CSX==NULL)
	{
		std::cerr << "CartOperator::Calc_EC: CSX not given or invalid!!!" << endl;
		return false;
	}

	MainOp->SetPos(0,0,0);

	// Same disjoint-X-slice parallel sampling as Operator_Multithread; CSXCAD
	// queries are read-only apart from its idempotent primitive used-flag.
	const unsigned int workers = std::min(GetSetupThreads(), numLines[0]);
	ParallelRanges(numLines[0], workers, [this](unsigned int start, unsigned int stop) {
		Calc_EC_Range(start, stop-1);
	});
	cout << "Metal: material EC threads: " << workers << endl;
	return true;
}

void Operator_Metal::CalcOperatorCoefficients()
{
	// Arithmetic only: reads EC arrays and writes disjoint X slabs. Each worker
	// uses its own AdrOp, since Calc_ECOperatorPos mutates the shared MainOp.
	const unsigned int workers = std::min(GetSetupThreads(), numLines[0]);
	ParallelRanges(numLines[0], workers, [this](unsigned int start, unsigned int stop) {
		AdrOp address(MainOp);
		unsigned int pos[3];
		for (pos[0]=start; pos[0]<stop; ++pos[0])
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2])
				{
					const unsigned int index = address.SetPos(pos[0],pos[1],pos[2]);
					for (int n=0; n<3; ++n) Calc_ECOperatorIndex(n,pos,index);
				}
	});
	cout << "Metal: operator coefficient threads: " << workers << endl;
}

bool Operator_Metal::CanReleaseECBeforeExtensions() const
{
	// Exact types, not derived types: future extensions must opt into this audit.
	// In particular series RLC, dispersive and conducting-sheet extensions need EC.
	for (const auto* extension : m_Op_exts)
		if (typeid(*extension) != typeid(Operator_Ext_UPML) &&
		    typeid(*extension) != typeid(Operator_Ext_Excitation) &&
		    typeid(*extension) != typeid(Operator_Ext_Mur_ABC))
		{
			cout << "Metal: keeping EC arrays resident for extension '" << extension->GetExtensionName() << "'" << endl;
			return false;
		}
	return true;
}

Engine* Operator_Metal::CreateEngine()
{
	m_Engine = Engine_Metal::New(this);
	return m_Engine;
}
