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
#include <cstdlib>
#include <typeinfo>
#include <thread>
#include <exception>
#include <algorithm>

using std::cout;
using std::endl;

Operator_Metal* Operator_Metal::New(unsigned int threads)
{
	cout << "Create FDTD operator (Metal field updates)" << endl;
	Operator_Metal* op = new Operator_Metal();
	op->Init();
	op->m_setupThreads = threads;
	return op;
}

Operator_Metal::Operator_Metal() : Operator_sse(), m_setupThreads(0)
{
}

bool Operator_Metal::Calc_EC()
{
	if (CSX==NULL)
	{
		std::cerr << "CartOperator::Calc_EC: CSX not given or invalid!!!" << endl;
		return false;
	}

	MainOp->SetPos(0,0,0);

	// The CPU engine parallelizes the material/geometry sampling by X range in
	// Operator_Multithread. Operator_Metal derives from Operator_sse, so without
	// this override the same CSXCAD-bound sampling runs on a single thread and
	// dominates Metal operator setup. Writes are disjoint per X slice and the
	// CSXCAD queries are read-only, matching the audited multithreaded path.
	// (CSXCAD still writes its idempotent used-flag on the winning primitive,
	// exactly as the existing Operator_Multithread path does.)
	unsigned int workers = GetSetupThreads() ? GetSetupThreads() : std::max(1U,std::thread::hardware_concurrency());
	workers = std::min(workers,numLines[0]);
	std::vector<std::thread> threads;
	std::vector<std::exception_ptr> errors(workers);
	auto run = [&](unsigned int worker) {
		try {
			const unsigned int start = numLines[0]*worker/workers;
			const unsigned int stop = numLines[0]*(worker+1)/workers;
			if (start < stop) Calc_EC_Range(start, stop-1);
		} catch (...) { errors[worker]=std::current_exception(); }
	};
	try {
		for (unsigned int i=0;i<workers;++i) threads.emplace_back(run,i);
	} catch (...) {
		for (auto& thread:threads) thread.join();
		throw;
	}
	for (auto& thread:threads) thread.join();
	for (auto error:errors) if (error) std::rethrow_exception(error);

	cout << "Metal: material EC threads: " << workers << endl;
	return true;
}

void Operator_Metal::CalcOperatorCoefficients()
{
	const char* serial = std::getenv("OPENEMS_METAL_SERIAL_COEFFICIENTS");
	if (serial && serial[0]=='1') { Operator::CalcOperatorCoefficients(); return; }
	unsigned int workers = GetSetupThreads() ? GetSetupThreads() : std::max(1U,std::thread::hardware_concurrency());
	workers = std::min(workers,numLines[0]);
	// This arithmetic pass touches no CSXCAD geometry: it only reads EC arrays
	// and writes disjoint X slabs of the operator coefficients. A thread-local
	// AdrOp computes the linear index, so the shared MainOp position is never
	// mutated concurrently.
	std::vector<std::thread> threads;
	std::vector<std::exception_ptr> errors(workers);
	auto run = [&](unsigned int worker) {
		try {
			AdrOp address(MainOp);
			unsigned int pos[3];
			for (pos[0]=numLines[0]*worker/workers; pos[0]<numLines[0]*(worker+1)/workers; ++pos[0])
				for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
					for (pos[2]=0; pos[2]<numLines[2]; ++pos[2]) {
						unsigned int index=address.SetPos(pos[0],pos[1],pos[2]);
						// Call the index form directly: Calc_ECOperatorPos would use the
						// shared MainOp position and is therefore not thread-safe.
						for (int n=0; n<3; ++n) Calc_ECOperatorIndex(n,pos,index);
					}
		} catch (...) { errors[worker]=std::current_exception(); }
	};
	try {
		for (unsigned int i=0;i<workers;++i) threads.emplace_back(run,i);
	} catch (...) {
		for (auto& thread:threads) thread.join();
		throw;
	}
	for (auto& thread:threads) thread.join();
	for (auto error:errors) if (error) std::rethrow_exception(error);
	cout << "Metal: operator coefficient threads: " << workers << endl;
}

bool Operator_Metal::CanReleaseECBeforeExtensions() const
{
	const char* setting = std::getenv("OPENEMS_METAL_EARLY_EC_FREE");
	if (setting && setting[0] == '0') return false;
	// Exact types, not derived types: future extensions must opt into this audit.
	// In particular series RLC, dispersive and conducting-sheet extensions need EC.
	for (const auto* extension : m_Op_exts)
		if (typeid(*extension) != typeid(Operator_Ext_UPML) &&
		    typeid(*extension) != typeid(Operator_Ext_Excitation) &&
		    typeid(*extension) != typeid(Operator_Ext_Mur_ABC)) return false;
	return true;
}

Engine* Operator_Metal::CreateEngine()
{
	m_Engine = Engine_Metal::New(this);
	return m_Engine;
}
