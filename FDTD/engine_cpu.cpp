/*
*	Copyright (C) 2026 Thorsten Liebig (Thorsten.Liebig@gmx.de)
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
#include "engine_cpu.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <typeinfo>
#ifdef OPENEMS_CPU_AVX2
#if defined(_MSC_VER) && !defined(__clang__)
#include <intrin.h>   // __cpuid(), __cpuidex(), _xgetbv()
#else
#include <cpuid.h>    // GCC and Clang, clang-cl too
#endif
#endif
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif
#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

#include "tools/global.h"
#include "tools/useful.h"
#include "extensions/engine_extension.h"
#include "extensions/engine_ext_upml.h"
#include "extensions/engine_ext_lorentzmaterial.h"
#include "extensions/engine_ext_debyematerial.h"
#include "extensions/engine_ext_absorbing_bc.h"
#include "extensions/operator_ext_absorbing_bc.h"
#include "extensions/operator_ext_cylinder.h"

using std::cout;
using std::endl;

static const char AXIS_LETTER[3] = {'I', 'J', 'K'};

Engine_CPU* Engine_CPU::New(const Operator* op, unsigned int numThreads)
{
	cout << "Create FDTD engine (" << GetTypeName(CPU) << ")" << endl;
	Engine_CPU* e = new Engine_CPU(op);
	e->SetNumThreads(numThreads);
	e->Init();
	return e;
}

Engine_CPU::Engine_CPU(const Operator* op) : Engine(op)
{
	m_type = CPU;
	m_FastAxis = 2;
	m_SplitAxis = 0;
	m_MidAxis = 1;
	m_Layout = "NIJK";
}

Engine_CPU::~Engine_CPU()
{
	Reset();
}

void Engine_CPU::Init()
{
	// OPENEMS_CPU_INERT=0: update the inert currents too
	const char* inert = std::getenv("OPENEMS_CPU_INERT");
	m_Model.reset(new CellModel(Op));
	m_Model->Analyse(!(inert && (std::string(inert)=="0")));
	ChooseLayout();
	Engine::Init();   // the fields in m_Layout (GetFieldLayout()), the extensions
	DropFusedExtensions();
	m_Model->Build(m_FastAxis, m_SplitAxis, m_MidAxis);
	if (Op->GetDebugFlags() & Operator::debugOperator)
		m_Model->Dump("cellmodel_dump");
	m_Model->FreeFused();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			m_State[f][n].assign(m_Model->Component(f, n).state, 0);
			const std::vector<double>& g = m_Model->G(f, n, m_FastAxis);
			m_Gf[f][n].assign(g.begin(), g.end());
		}
	InitKernels();
}

void Engine_CPU::DropFusedExtensions()
{
	// the engine updates their cells itself (the classes UPML, UPML_PROFILE, DISP)
	std::vector<Engine_Extension*> keep;
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
	{
		Engine_Extension* e = m_Eng_exts[n];
		if (!m_Model->IsFused(e->GetOperatorExtension()))
		{
			keep.push_back(e);
			continue;
		}
		const std::type_info& type = typeid(*e);
		if ((type!=typeid(Engine_Ext_UPML)) && (type!=typeid(Engine_Ext_LorentzMaterial)) && (type!=typeid(Engine_Ext_DebyeMaterial)))
			throw std::logic_error("Engine_CPU: an unknown engine extension of the fused " + e->GetExtensionName());
		delete e;
	}
	m_Eng_exts = keep;
}

void Engine_CPU::Reset()
{
	Engine::Reset();
	m_Model.reset();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			m_State[f][n].clear();
			m_Gf[f][n].clear();
		}
}

void Engine_CPU::ChooseLayout()
{
	m_FastAxis = m_Model->RecommendFastAxis();
	// the box plans: a closed alpha mesh couples the first and the last alpha plane (Engine_Ext_Cylinder),
	// a local absorber its sheet and the next line along its normal within a half step; as the fast
	// axis, which every box spans, the blocked plans can run them (BuildBlockedPlans())
	bool alpha = false;
	int normal = -1;
	for (size_t n=0; n<Op->GetNumberOfExtentions(); ++n)
	{
		const Operator_Extension* e = Op->GetExtension(n);
		const Operator_Ext_Cylinder* cyl = dynamic_cast<const Operator_Ext_Cylinder*>(e);
		if (cyl && cyl->ClosedAlpha())
			alpha = true;
		const Operator_Ext_Absorbing_BC* abc = dynamic_cast<const Operator_Ext_Absorbing_BC*>(e);
		if (abc && (normal<0))
			normal = abc->GetNormalAxis();
		else if (abc && (normal!=abc->GetNormalAxis()))
			normal = 3;   // several: no blocking
	}
	if (alpha)
		m_FastAxis = 1;
	else if ((normal>=0) && (normal<3))   // the cylinder's alpha first: its box hooks need it
		m_FastAxis = normal;
	const char* env = std::getenv("OPENEMS_CPU_FAST_AXIS");
	if (env && env[0])
	{
		const char c = env[0];
		if ((c>='x') && (c<='z'))
			m_FastAxis = c-'x';
		else
			throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_FAST_AXIS must be x, y or z");
	}

	// the threads split the slow axis with more lines
	const int a1 = (m_FastAxis+1)%3, a2 = (m_FastAxis+2)%3;
	m_SplitAxis = (numLines[a1]>=numLines[a2]) ? a1 : a2;
	m_MidAxis = 3 - m_FastAxis - m_SplitAxis;
	m_Layout = std::string("N") + AXIS_LETTER[m_SplitAxis] + AXIS_LETTER[m_MidAxis] + AXIS_LETTER[m_FastAxis];
	if (g_settings.GetVerboseLevel()>0)
		cout << "Engine_CPU: fast axis " << (char)('x'+m_FastAxis) << ", layout " << m_Layout << endl;
}

// ----------------------------------------------------------------------------
// the kernels are in engine_cpu_kernels.inc, compiled per instruction set

#ifdef OPENEMS_CPU_AVX2
// AVX2 and FMA, and the AVX registers saved by the system: from CPUID and XGETBV directly.
// __builtin_cpu_supports() needs __cpu_model of the compiler runtime, which e.g. clang
// on Windows does not link (x86-64 only, see CMakeLists.txt)
static bool CpuHasAVX2FMA()
{
	unsigned int b, c;   // CPUID leaf 1 ECX, leaf 7 EBX
	unsigned long long xcr0;
#if defined(_MSC_VER) && !defined(__clang__)
	int r[4];
	__cpuid(r, 0);
	if (r[0]<7)
		return false;
	__cpuid(r, 1);
	c = (unsigned int)r[2];
	if (!((c>>27) & 1))   // OSXSAVE: XGETBV available
		return false;
	xcr0 = _xgetbv(0);
	__cpuidex(r, 7, 0);
	b = (unsigned int)r[1];
#else
	unsigned int a, d, c7;
	if (!__get_cpuid(1, &a, &b, &c, &d))
		return false;
	if (!((c>>27) & 1))   // OSXSAVE: XGETBV available
		return false;
	unsigned int lo, hi;
	__asm__ volatile ("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
	xcr0 = ((unsigned long long)hi<<32) | lo;
	if (!__get_cpuid_count(7, 0, &a, &b, &c7, &d))   // c keeps leaf 1
		return false;
#endif
	const bool fma = (c>>12) & 1, avx = (c>>28) & 1;
	return fma && avx && ((xcr0 & 6)==6)   // the SSE and the AVX state saved
	       && ((b>>5) & 1);                // AVX2
}
#endif

void Engine_CPU::InitKernels()
{
	Engine_CPU_Kernel& k = m_Kernel;
	k.F = m_FastAxis;
	k.S = m_SplitAxis;
	k.M = m_MidAxis;
	for (int a=0; a<3; ++a)
	{
		k.numLines[a] = numLines[a];
		k.s[a] = volt_ptr->stride(1+a);
	}
	k.sN = volt_ptr->stride(0);
	k.field[0] = volt_ptr->data();
	k.field[1] = curr_ptr->data();
	for (int f=0; f<2; ++f)
		for (int n=0; n<3; ++n)
		{
			Engine_CPU_Kernel::Component& c = k.comp[f][n];
			const CellModel_Component& mc = m_Model->Component(f, n);
			c.runs = mc.runs.data();
			c.first = mc.first.data();
			c.sets = mc.sets.data();
			c.state = m_State[f][n].data();
			c.poles = mc.poles.data();
			for (int q=0; q<5; ++q)
				c.profile[q] = mc.profile[q].data();
			c.gf = m_Gf[f][n].data();
			c.gS = m_Model->G(f, n, m_SplitAxis).data();
			c.gM = m_Model->G(f, n, m_MidAxis).data();
		}

	// the instruction set: AVX2 and FMA if the CPU has both (OPENEMS_CPU_ISA=generic overrides)
	m_ISA = "generic";
	m_UpdateVoltages = Engine_CPU_Generic::UpdateVoltages;
	m_UpdateCurrents = Engine_CPU_Generic::UpdateCurrents;
	const char* env = std::getenv("OPENEMS_CPU_ISA");
	const std::string want = env ? env : "";
	if ((want!="") && (want!="generic") && (want!="avx2"))
		throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_ISA must be generic or avx2");
#ifdef OPENEMS_CPU_AVX2
	if (want!="generic")
	{
		if (CpuHasAVX2FMA())
		{
			m_ISA = "avx2";
			m_UpdateVoltages = Engine_CPU_AVX2::UpdateVoltages;
			m_UpdateCurrents = Engine_CPU_AVX2::UpdateCurrents;
		}
		else if (want=="avx2")
			throw std::runtime_error("Engine_CPU: OPENEMS_CPU_ISA=avx2, but the CPU lacks AVX2 or FMA");
	}
#else
	if (want=="avx2")
		throw std::runtime_error("Engine_CPU: OPENEMS_CPU_ISA=avx2, but this build has no AVX2 kernels");
#endif
	if (g_settings.GetVerboseLevel()>0)
		cout << "Engine_CPU: kernels for " << m_ISA << endl;
}

// the bytes of the last level cache the cores share, and where they come from; 32 MiB if the system does not say
static size_t CacheBytes(std::string& source)
{
#if defined(_WIN32)
	// the largest cache of the highest level
	DWORD len = 0;
	GetLogicalProcessorInformationEx(RelationCache, NULL, &len);
	if (len>0)
	{
		std::vector<char> buf(len);
		SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* info = (SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)buf.data();
		if (GetLogicalProcessorInformationEx(RelationCache, info, &len))
		{
			int level = 0;
			size_t bytes = 0;
			for (DWORD off=0; off<len; )
			{
				const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX* e = (const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*)(buf.data()+off);
				if ((e->Relationship==RelationCache) && (e->Cache.Type!=CacheInstruction))
				{
					const int l = e->Cache.Level;
					if ((l>level) || ((l==level) && (e->Cache.CacheSize>bytes)))
					{
						level = l;
						bytes = e->Cache.CacheSize;
					}
				}
				off += e->Size;
			}
			if (bytes>0)
			{
				source = "L" + std::to_string(level) + " (Windows)";
				return bytes;
			}
		}
	}
#elif defined(__APPLE__)
	// an Intel Mac has an L3; Apple silicon has none (its system cache is not the cores'): the L2 of the performance cores
	const char* names[] = {"hw.l3cachesize", "hw.perflevel0.l2cachesize", "hw.l2cachesize"};
	for (const char* name : names)
	{
		int64_t bytes = 0;
		size_t size = sizeof(bytes);
		if ((sysctlbyname(name, &bytes, &size, NULL, 0)==0) && (bytes>0))
		{
			source = name;
			return (size_t)bytes;
		}
	}
#else
#ifdef _SC_LEVEL3_CACHE_SIZE
	const long l3 = sysconf(_SC_LEVEL3_CACHE_SIZE);
	if (l3>0)
	{
		source = "L3 (sysconf)";
		return l3;
	}
#endif
	// some systems (e.g. ARM) do not say through sysconf: the cache descriptions of cpu0, the highest level
	size_t best = 0;
	int bestLevel = 0;
	for (int index=0; index<8; ++index)
	{
		const std::string dir = "/sys/devices/system/cpu/cpu0/cache/index" + std::to_string(index) + "/";
		std::ifstream fl(dir + "level"), ft(dir + "type"), fs(dir + "size");
		int level = 0;
		std::string type, size;
		if (!(fl >> level) || !(ft >> type) || !(fs >> size) || (type=="Instruction") || size.empty())
			continue;
		size_t bytes = std::strtoul(size.c_str(), NULL, 10);
		const char unit = size.back();
		if ((unit=='K') || (unit=='k'))
			bytes <<= 10;
		else if ((unit=='M') || (unit=='m'))
			bytes <<= 20;
		if ((level>bestLevel) && (bytes>0))
		{
			best = bytes;
			bestLevel = level;
		}
	}
	if (best>0)
	{
		source = "L" + std::to_string(bestLevel) + " (sysfs)";
		return best;
	}
#endif
	source = "assumed, the system does not say";
	return 32u<<20;
}

bool Engine_CPU::BuildBlockedPlans(unsigned int numThreads, std::vector<Engine_BlockedPlan>& plans)
{
	const char* env = std::getenv("OPENEMS_CPU_TBLOCK");
	std::string why;
	unsigned int k = 0, W = 0;
	if (env && (std::string(env)=="0"))
		why = "OPENEMS_CPU_TBLOCK=0";
	else if (m_SubGrid)
		why = "a cylindrical multi-grid";
	for (size_t n=0; n<m_Eng_exts.size() && why.empty(); ++n)
	{
		const Engine_Ext_Absorbing_BC* abc = dynamic_cast<const Engine_Ext_Absorbing_BC*>(m_Eng_exts[n]);
		if (!m_Eng_exts[n]->SupportsBoxes())
			why = "the " + m_Eng_exts[n]->GetExtensionName() + " has no box hooks";
		else if (abc && (abc->GetNormalAxis()!=m_FastAxis))
			why = "a local absorber not normal to the fast axis";
	}
	if (why.empty() && env && env[0])
	{
		char* end = NULL;
		k = (unsigned int)std::strtoul(env, &end, 10);
		if ((end!=NULL) && (*end==':'))
			W = (unsigned int)std::strtoul(end+1, &end, 10);
		if ((k==0) || (end==NULL) || (*end!='\0'))
			throw std::invalid_argument("Engine_CPU: OPENEMS_CPU_TBLOCK must be 0, k or k:W");
		if (W==0)
			W = 2*k+1;
	}
	else if (why.empty())
	{
		// the tile: the fields, the UPML and ADE state and the runs of its planes in the cache. The timesteps
		// per block from the whole cache, the width then 3/4 of it as long as it holds them (measured: a tile
		// that fills the whole cache is much slower, the rest of the working set evicts it; but fewer
		// timesteps per block cost more than that)
		size_t runs = 0, state = 0;
		for (int f=0; f<2; ++f)
			for (int n=0; n<3; ++n)
			{
				runs += m_Model->Component(f, n).runs.size();
				state += m_State[f][n].size();
			}
		double ext = 0;   // the per-cell state of the extensions (dispersive materials)
		for (size_t n=0; n<m_Eng_exts.size(); ++n)
			ext += m_Eng_exts[n]->BoxBytes();
		const double plane = (double)numLines[m_MidAxis]*numLines[m_FastAxis];
		const double fieldBytes = 6*sizeof(FDTD_FLOAT);   // the voltages and currents of a node
		const double bytes = plane*fieldBytes + (sizeof(float)*state + sizeof(CellModel_Run)*runs + ext)/numLines[m_SplitAxis];
		std::string source;
		const size_t cache = CacheBytes(source);
		if (g_settings.GetVerboseLevel()>0)
			cout << "Engine_CPU: cache " << cache/1048576.0 << " MiB, " << source << endl;
		const unsigned int full = (unsigned int)(cache/bytes);
		k = std::min(8u, (full>0) ? (full-1)/2 : 0);
		W = std::max(2*k+1, (unsigned int)(0.75*full));
		if (k<2)
			why = "a plane is too large for the cache";
		else if (numLines[m_SplitAxis] < 2*W)
			why = "the model fits the cache";
	}
	if (!why.empty())
	{
		if (g_settings.GetVerboseLevel()>0)
			cout << "Engine_CPU: no temporal blocking: " << why << endl;
		return false;
	}
	for (size_t n=0; n<m_Eng_exts.size(); ++n)
		m_Eng_exts[n]->PrepareBoxes(numLines);
	// every k down to 2: the processing between two calls of IterateTS() may leave fewer timesteps
	for (unsigned int kk=k; kk>=2; --kk)
		plans.push_back(Engine_BlockedPlan::Trapezoids(numLines, numThreads, m_SplitAxis, m_MidAxis, W, kk));
	if (g_settings.GetVerboseLevel()>0)
		cout << "Engine_CPU: temporal blocking, up to " << k << " timesteps per block, tiles of " << W << " lines along "
		     << (char)('x'+m_SplitAxis) << ", " << plans[0].phases.size() << " phases" << endl;
	return true;
}

void Engine_CPU::UpdateVoltages(const Engine_Box& box)
{
	m_UpdateVoltages(m_Kernel, box.start, box.stop);
}

void Engine_CPU::UpdateCurrents(const Engine_Box& box)
{
	m_UpdateCurrents(m_Kernel, box.start, box.stop);
}
