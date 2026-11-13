/*
*	Copyright (C) 2026 Sean Mollet (sean@malmoset.com)
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

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <stdexcept>

#include "engine_gpu.h"
#include "cell_model.h"
#include "gpu_cell_model.h"
#include "gpu_ext.h"
#include "gpu_backend.h"
#include "engine_subgrid.h"
#include "extensions/engine_extension.h"
#include "extensions/engine_ext_absorbing_bc.h"
#include "extensions/operator_extension.h"

using std::cout;
using std::cerr;
using std::endl;

//! \brief construct an Engine_GPU instance
//! it's the responsibility of the caller to free the returned pointer
Engine_GPU* Engine_GPU::New(const Operator* op, GPU_Backend* backend, GPU_Backend* parent, const std::function<void()>& modelBuilt)
{
	Engine_GPU* e = new Engine_GPU(op, backend, parent);
	e->m_ModelBuilt = modelBuilt;
	try
	{
		e->Init();
	}
	catch (...)
	{
		delete e;
		throw;
	}
	return e;
}

Engine_GPU::Engine_GPU(const Operator* op, GPU_Backend* backend, GPU_Backend* parent) : Engine(op)
{
	m_ParentBackend = parent;
	m_type = GPU;
	// ours from here on, Reset() frees it even if Init() fails
	m_Backend = backend;
	m_SharedMemory = false;
	m_SnapshotsUsed = false;
	m_SnapshotGatherSet = false;
	m_SnapshotGatherSize[0] = m_SnapshotGatherSize[1] = 0;
	m_FieldDFTUsed = false;
	m_DumpPath = 0;
	m_StaleVolt.stale = m_StaleCurr.stale = false;
	m_StaleVolt.full_this_batch = m_StaleCurr.full_this_batch = false;
	m_StaleVolt.full_last_batch = m_StaleCurr.full_last_batch = false;
}

Engine_GPU::~Engine_GPU()
{
	Reset();
}

void Engine_GPU::Init()
{
	// a sub-grid shares the device work stream of its parent grid
	if (m_ParentBackend)
		m_Backend = m_ParentBackend->NewSubGridBackend();
	if (m_Backend==NULL)
		throw std::runtime_error("Engine_GPU::Init: no GPU backend");
	cout << "Create FDTD engine (GPU, backend: " << m_Backend->GetName() << ")" << endl;
	m_DumpPath = 0;
	if (const char* env = getenv("OPENEMS_GPU_DUMPS"))
	{
		const std::string path(env);
		m_DumpPath = (path=="snapshot") ? 1 : ((path=="host") ? 2 : 0);
		if (m_DumpPath)
			cout << "Engine_GPU: OPENEMS_GPU_DUMPS=" << path << ": field dumps " << (m_DumpPath==1 ? "without device sums" : "on the host") << endl;
	}
	// the cell model of the CPU engine along z (the ArrayNIJK layout); its fused extensions
	// (UPML, dispersive materials) are part of the main updates, their engine extensions go.
	// The memory peak as in Engine_CPU::Init(): the operator's coefficients go once the model
	// is built, the model is expanded and uploaded (and freed) before the host mirror comes
	CellModel model(Op);
	model.Analyse(true);
	model.Build(2, 0, 1);
	if (Op->GetDebugFlags() & Operator::debugOperator)
		model.Dump("cellmodel_dump");
	if (m_ModelBuilt)
	{
		m_ModelBuilt();
		m_ModelBuilt = nullptr;
	}
	{
		const GPU_CellModel expanded(model, numLines);
		m_NodeBytes = expanded.NodeBytes();
		if (!m_Backend->Init(Op, expanded))
			throw std::runtime_error("Engine_GPU::Init: GPU backend initialization failed");
	}
	// allocates the host mirror and creates the extensions
	Engine::Init();
	model.DropFused(m_Eng_exts);
	model.FreeFused();

	// with unified memory the host mirror is a view of the device fields
	m_SharedMemory = m_Backend->GetShared(false) && m_Backend->GetShared(true);
	if (m_SharedMemory)
	{
		volt_ptr->InitView("volt", {numLines[0], numLines[1], numLines[2]}, m_Backend->GetShared(false));
		curr_ptr->InitView("curr", {numLines[0], numLines[1], numLines[2]}, m_Backend->GetShared(true));
	}

	// every extension on the device (a device implementation of each, gpu_ext.h): the fields stay there
	std::string missing;
	for (Engine_Extension* eng_ext : m_Eng_exts)
	{
		GPU_Extension* gpu_ext = GPU_CreateExtension(eng_ext, m_Backend, numLines);
		if (gpu_ext)
			m_GPU_exts.push_back(gpu_ext);
		else
			missing += (missing.empty() ? "" : ", ") + eng_ext->GetExtensionName();
	}
	if (!missing.empty())
		throw std::runtime_error("Engine_GPU: no device implementation (gpu_ext.h) of the " + missing);
	m_BlockingPrepared = false;
}

void Engine_GPU::PrepareBlocking()
{
	m_BlockCandidates.clear();
	m_Block = -1;
	m_Trial = -1;
	m_FlatSplit = false;
	// as Engine_CPU: tiles given (no trials then, unless asked for), or the cache size to choose them from
	std::vector< std::pair<unsigned int, unsigned int> > tiles;   // (k, W)
	const bool explicitTile = Engine_TilesFromEnv("OPENEMS_GPU_TBLOCK", tiles);
	m_TrialSeconds = Engine_NumberFromEnv("OPENEMS_GPU_TRIAL", tiles.empty() ? 0.3 : 0, true);
	std::string why;
	if (explicitTile && tiles.empty())
		why = "OPENEMS_GPU_TBLOCK=0";
	else if (m_SubGrid || m_ParentBackend)
		why = "a cylindrical multi-grid";
	// a local absorber needs its sheet and the next line along its normal in one box: the boxes split x
	// (Engine_CPU blocks along its normal instead, Engine_CPU::ChooseLayout())
	for (size_t n=0; n<m_Eng_exts.size() && why.empty(); ++n)
	{
		const Engine_Ext_Absorbing_BC* abc = dynamic_cast<const Engine_Ext_Absorbing_BC*>(m_Eng_exts[n]);
		if (abc && (abc->GetNormalAxis()==0))
			why = "a local absorber normal to x";
	}

	std::string skipFused;   // why the fused step is not worth trying
	if (why.empty() && !explicitTile)
	{
		// the tile: what the main updates read of its planes (GPU_CellModel::NodeBytes()) in the last-level cache. Of its
		// size, if the backend tells it (or OPENEMS_GPU_CACHE): tiles of 1/2, 3/4 and all of it; else of 8 to
		// 256 MiB. The trials keep the fastest
		std::string source;
		double cache = (double)m_Backend->CacheBytes(source);
		if (const double mib = Engine_NumberFromEnv("OPENEMS_GPU_CACHE", 0, false))
		{
			cache = mib*1048576.0;
			source = "OPENEMS_GPU_CACHE";
		}
		std::vector<double> budgets;
		if (cache>0)
			budgets = {0.5*cache, 0.75*cache, cache};
		else
			for (double mib=8; mib<=256; mib*=2)
				budgets.push_back(mib*1048576.0);
		if (g_settings.GetVerboseLevel()>0)
		{
			if (cache>0)
				cout << "Engine_GPU: last-level cache " << cache/1048576.0 << " MiB, " << source << endl;
			else
				cout << "Engine_GPU: last-level cache unknown, tiles for 8 to 256 MiB" << endl;
		}
		const double plane = (double)numLines[1]*numLines[2]*m_NodeBytes;
		for (double bytes : budgets)
		{
			const unsigned int W = (unsigned int)(bytes/plane);
			const unsigned int k = std::min(8u, (W>0) ? (W-1)/2 : 0);
			if ((k>=2) && (2*W<=numLines[0]) && (tiles.empty() || (tiles.back().second!=W)))
				tiles.push_back(std::make_pair(k, W));
		}
		if (tiles.empty())
			why = "the model fits the cache, or a plane is too large for it";
		// a tile of the full depth streams the fields once per 8 timesteps, the fused step once per timestep:
		// with the cache known to hold one, the fused step cannot win, so neither its memory nor its trial
		else if ((cache>0) && (tiles.back().first>=8))
			skipFused = "the cache holds tiles of 8 timesteps";
	}
	PrepareFused(skipFused);
	// the flat sweep: the fused step against the two updates, unless OPENEMS_GPU_FUSED=1 asks for it
	m_FlatTrial = m_Fused && !m_FusedForced;
	if (!why.empty())
	{
		if (g_settings.GetVerboseLevel()>0)
			cout << "Engine_GPU: no temporal blocking: " << why << endl;
		if (m_FlatTrial && (m_TrialSeconds>0))
		{
			if (g_settings.GetVerboseLevel()>0)
				cout << "Engine_GPU: trying the fused step and the two updates" << endl;
			StartTrials();
		}
		return;
	}
	for (const std::pair<unsigned int, unsigned int>& t : tiles)
	{
		std::vector<Engine_BlockedPlan> plans;
		for (unsigned int k=t.first; k>=2; --k)
			plans.push_back(Engine_BlockedPlan::Parallelograms(numLines, 1, 0, 1, t.second, k));
		m_BlockCandidates.push_back(plans);
	}
	if (g_settings.GetVerboseLevel()>0)
	{
		cout << "Engine_GPU: temporal blocking, tiles to try (timesteps:lines)";
		for (size_t c=0; c<tiles.size(); ++c)
			cout << " " << c+1 << ") " << tiles[c].first << ":" << tiles[c].second << " (" << m_BlockCandidates[c].front().tiles << " tiles)";
		cout << (m_TrialSeconds>0 ? (m_FlatTrial ? ", the fused step and the two updates" : ", and the flat sweep") : "") << endl;
	}
	m_Block = 0;
	if ((m_TrialSeconds<=0) && m_Fused)
	{
		// a tile without trials: the flat steps are rare
		m_Backend->ReleaseFused();
		m_Fused = false;
	}
	if (m_TrialSeconds>0)
		StartTrials();
}

void Engine_GPU::StartTrials()
{
	m_Trial = 0;
	m_TrialSteps = 0;
	m_TrialSpeed.assign(m_BlockCandidates.size() + 1 + (m_FlatTrial ? 1 : 0), 0);
	m_TrialTime = 0;
	m_TrialMark = -1;
	m_TrialWarm = true;
	m_Block = m_BlockCandidates.empty() ? -1 : 0;
	m_FlatSplit = false;
}

void Engine_GPU::PrepareFused(const std::string& skip)
{
	m_Fused = false;
	const char* env = getenv("OPENEMS_GPU_FUSED");
	m_FusedForced = env && (std::string(env)=="1");
	if ((env && (std::string(env)=="0")) || m_SubGrid || m_ParentBackend)
		return;
	if (!skip.empty() && !m_FusedForced)
	{
		if (g_settings.GetVerboseLevel()>0)
			cout << "Engine_GPU: no fused step: " << skip << endl;
		return;
	}
	// the currents that wait for the hooks between the main updates: those the hooks read or
	// write, and those whose curl reads a voltage they write (the voltages at the node and
	// above: a current of node i reads the voltages at i, i+sx, i+sy and i+sz)
	std::vector<uint64_t> volt, late;
	for (const GPU_Extension* ext : m_GPU_exts)
		if (!ext->MidStepNodes(volt, late))
		{
			if (g_settings.GetVerboseLevel()>0)
				cout << "Engine_GPU: no fused step: an extension does not tell the nodes of its hooks" << endl;
			return;
		}
	const uint64_t sz = 1, sy = numLines[2], sx = (uint64_t)numLines[1]*numLines[2];
	for (uint64_t v : volt)
	{
		late.push_back(v);
		if (v/sx > 0)
			late.push_back(v-sx);
		if ((v/sy)%numLines[1] > 0)
			late.push_back(v-sy);
		if (v%numLines[2] > 0)
			late.push_back(v-sz);
	}
	std::sort(late.begin(), late.end());
	late.erase(std::unique(late.begin(), late.end()), late.end());
	m_Fused = m_Backend->PrepareFused(late);
	if (m_Fused && (g_settings.GetVerboseLevel()>0))
		cout << "Engine_GPU: fused step, " << late.size() << " late currents" << endl;
}

template <typename Update>
void Engine_GPU::DeviceHalfStep(bool currents, const Engine_Box& box, unsigned int ts, Update update)
{
	// the extensions see the last line of an axis too, which has no current update (as Engine_ThreadPool::BlockedStep())
	Engine_Box hook = box;
	if (currents)
		for (int a=0; a<3; ++a)
			if (hook.stop[a]+1==numLines[a])
				hook.stop[a] = numLines[a];
	for (size_t n=m_GPU_exts.size(); n-->0; )
		m_GPU_exts[n]->Hook(currents, GPU_PRE, hook, ts);
	update();
	for (GPU_Extension* ext : m_GPU_exts)
		ext->Hook(currents, GPU_POST, hook, ts);
	for (GPU_Extension* ext : m_GPU_exts)
		ext->Hook(currents, GPU_APPLY, hook, ts);
}

void Engine_GPU::FusedStep()
{
	DeviceHalfStep(false, Domain(), numTS, [this]() {m_Backend->FusedUpdate();});
	DeviceHalfStep(true, Domain(), numTS, [this]() {m_Backend->SeamUpdate();});
}

unsigned int Engine_GPU::BlockLength(unsigned int left, unsigned int kmax) const
{
	kmax = std::min(left, kmax);
	unsigned int k = 0;
	while (k<kmax)
	{
		bool flat = false;
		for (const GPU_Extension* ext : m_GPU_exts)
			flat = flat || ext->NeedsFlatStep(numTS+k);
		if (flat)
			break;
		++k;
	}
	return k;
}

void Engine_GPU::RunBlock(const Engine_BlockedPlan& plan)
{
	for (const std::vector< std::vector<Engine_Step> >& phase : plan.phases)
		for (const Engine_Step& step : phase.at(0))
			BlockedStep(step, numTS + step.t);
	for (unsigned int t=0; t<plan.k; ++t)
		NextTimestep();
}

void Engine_GPU::BlockedStep(const Engine_Step& step, unsigned int ts)
{
	if (step.volt)
		DeviceHalfStep(false, step.box, ts, [&]() {m_Backend->UpdateVoltages(step.box);});
	else
		DeviceHalfStep(true, step.box, ts, [&]() {m_Backend->UpdateCurrents(step.box);});
}

static double WallSeconds()
{
	return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

void Engine_GPU::TrialClock(bool run)
{
	if ((m_Trial<0) || (run==(m_TrialMark>=0)))
		return;
	m_Backend->Synchronize();
	const double now = WallSeconds();
	if (run)
		m_TrialMark = now;
	else
	{
		m_TrialTime += now - m_TrialMark;
		m_TrialMark = -1;
	}
}

void Engine_GPU::TrialProgress(unsigned int steps)
{
	if (m_Trial<0)
		return;
	if (m_TrialWarm)
	{
		// the first steps of a candidate are not timed (the warm-up of its kernels and caches)
		TrialClock(false);
		m_TrialTime = 0;
		m_TrialSteps = 0;
		m_TrialWarm = false;
		TrialClock(true);
		return;
	}
	m_TrialSteps += steps;
	if (m_TrialTime + WallSeconds()-m_TrialMark < m_TrialSeconds)
		return;
	TrialClock(false);
	m_TrialSpeed[m_Trial] = m_TrialSteps/std::max(m_TrialTime, 1e-9);
	++m_Trial;
	m_TrialTime = 0;
	m_TrialWarm = true;
	// the candidates: the tiles, the flat sweep (the fused step, if prepared), the two updates (m_FlatTrial)
	const int tiles = m_BlockCandidates.size();
	if (m_Trial<(int)m_TrialSpeed.size())
	{
		TrialClock(true);
		m_Block = (m_Trial<tiles) ? m_Trial : -1;
		m_FlatSplit = (m_Trial==tiles+1);
		return;
	}
	// all tried: the fastest
	int best = 0;
	for (size_t c=1; c<m_TrialSpeed.size(); ++c)
		if (m_TrialSpeed[c]>m_TrialSpeed[best])
			best = c;
	m_Block = (best<tiles) ? best : -1;
	m_Trial = -1;
	m_FlatSplit = false;
	if (m_Fused && ((m_Block>=0) || (best==tiles+1)))
	{
		// the memory of the fused step back (the flat steps of the blocked runs are rare)
		m_Backend->ReleaseFused();
		m_Fused = false;
	}
	if (g_settings.GetVerboseLevel()>0)
	{
		const double cells = (double)numLines[0]*numLines[1]*numLines[2]/1e6;
		auto name = [&](int c) {
			if (c<tiles)
				return "tile " + std::to_string(c+1);
			if (!m_FlatTrial)
				return std::string("the flat sweep");
			return std::string((c==tiles) ? "the fused step" : "the two updates");
		};
		cout << "Engine_GPU: tried";
		for (size_t c=0; c<m_TrialSpeed.size(); ++c)
			cout << (c ? ", " : " ") << name(c) << " " << (int)(m_TrialSpeed[c]*cells) << " MC/s";
		cout << "; running " << name(best) << endl;
	}
}

void Engine_GPU::ClearGPUExtensions()
{
	for (size_t n=0; n<m_GPU_exts.size(); ++n)
		delete m_GPU_exts.at(n);
	m_GPU_exts.clear();
}

void Engine_GPU::Reset()
{
	// the sub-grid backend shares the work stream of ours, the coupling refers to both
	ClearSubGrid();
	// the device extensions refer to the backend and to the engine extensions
	ClearGPUExtensions();
	delete m_Backend;
	m_Backend = NULL;
	Engine::Reset();
}

void Engine_GPU::VoltagesToHost()
{
	if (m_SharedMemory)
		m_Backend->Synchronize();
	else
		m_Backend->Download(false, *volt_ptr);
	m_StaleVolt.stale = false;
	m_StaleVolt.lines.clear();
}

void Engine_GPU::CurrentsToHost()
{
	if (m_SharedMemory)
		m_Backend->Synchronize();
	else
		m_Backend->Download(true, *curr_ptr);
	m_StaleCurr.stale = false;
	m_StaleCurr.lines.clear();
}

FDTD_FLOAT Engine_GPU::GetVolt(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
{
	if (m_StaleVolt.stale)
		return ReadStale(false, n, x, y, z);
	return Engine::GetVolt(n, x, y, z);
}

FDTD_FLOAT Engine_GPU::GetCurr(unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
{
	if (m_StaleCurr.stale)
		return ReadStale(true, n, x, y, z);
	return Engine::GetCurr(n, x, y, z);
}

// Reading a z-line takes about the latency of a transfer, the whole field (24 bytes
// per node) a lot longer. Probes read a few lines per batch, dumps all of them.
static const size_t MAX_STALE_LINES = 256;

FDTD_FLOAT Engine_GPU::ReadStale(bool currents, unsigned int n, unsigned int x, unsigned int y, unsigned int z) const
{
	Engine_GPU* self = const_cast<Engine_GPU*>(this);   // the host mirror is a cache of the device fields
	StaleField& field = currents ? m_StaleCurr : m_StaleVolt;
	const size_t key = ((size_t)n*numLines[0] + x)*numLines[1] + y;
	std::unordered_map<size_t, std::vector<FDTD_FLOAT>>::const_iterator it = field.lines.find(key);
	if (it!=field.lines.end())
		return it->second.at(z);

	if (!field.full_last_batch && (field.lines.size()<MAX_STALE_LINES))
	{
		std::vector<FDTD_FLOAT> line(numLines[2]);
		m_Backend->DeviceDownload(line.data(), m_Backend->DeviceField(currents) + key*numLines[2], numLines[2]*sizeof(FDTD_FLOAT));
		const FDTD_FLOAT value = line.at(z);
		field.lines[key].swap(line);
		return value;
	}

	// many lines: the whole field
	if (currents)
		self->CurrentsToHost();
	else
		self->VoltagesToHost();
	field.full_this_batch = true;
	return currents ? Engine::GetCurr(n, x, y, z) : Engine::GetVolt(n, x, y, z);
}

void Engine_GPU::MarkStale(StaleField& field)
{
	// a batch that read the whole field is likely followed by another one (e.g. dumps)
	field.full_last_batch = field.full_this_batch;
	field.full_this_batch = false;
	field.stale = true;
	field.lines.clear();
}

bool Engine_GPU::SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (m_DumpPath>=2)
		return false;
	// not for a sub-grid: the dumps read the base grid, which FinishBatch() gave the
	// sub-grid fields on the work stream before
	if (m_ParentBackend)
		return false;
	// A snapshot has a fixed cost (a separate command buffer or stream, the background thread), which
	// small grids with short batches between dumps do not recover (e.g. 69k cells dumped
	// every few timesteps: twice the engine time). The copy itself scales with the grid,
	// like the timesteps between the dumps.
	if ((size_t)numLines[0]*numLines[1]*numLines[2] < 512*1024)
		return false;
	if (CanSnapshotGather() && !m_SnapshotGatherSet)
	{
		// once: without enough device memory the backend has no snapshots
		m_Backend->SetSnapshotGather(m_SnapshotGather[0], m_SnapshotGather[1]);
		m_SnapshotGatherSet = true;
		std::vector<Engine_GatherRegion>().swap(m_SnapshotGather[0]);
		std::vector<Engine_GatherRegion>().swap(m_SnapshotGather[1]);
	}
	if (!m_Backend->SnapshotFields(volt, curr))
		return false;
	if (!m_SnapshotsUsed && (g_settings.GetVerboseLevel()>1))
		cout << "Engine_GPU: field dumps from snapshots of the fields" << endl;
	m_SnapshotsUsed = true;
	return true;
}

void Engine_GPU::WaitSnapshot() const
{
	m_Backend->WaitSnapshot();
}

bool Engine_GPU::CanSnapshotGather() const
{
	if (m_ParentBackend || (m_DumpPath>=2))
		return false;
	return m_Backend->CanSnapshotGather();
}

bool Engine_GPU::AddSnapshotGather(bool h_field, const Engine_GatherRegion& region, size_t& offset)
{
	if (m_SnapshotGatherSet || !CanSnapshotGather())
		return false;
	const int f = h_field ? 1 : 0;
	offset = m_SnapshotGatherSize[f];
	m_SnapshotGather[f].push_back(region);
	m_SnapshotGatherSize[f] += region.Size();
	return true;
}

int Engine_GPU::AddFieldDFT(bool h_field, const Engine_GatherRegion& region, unsigned int count)
{
	if (m_DumpPath>=1)
		return -1;
	const int id = m_Backend->AddFieldDFT(h_field, region, count);
	if ((id>=0) && !m_FieldDFTUsed)
	{
		if (g_settings.GetVerboseLevel()>1)
			cout << "Engine_GPU: frequency domain dumps accumulated on the device" << endl;
		m_FieldDFTUsed = true;
	}
	return id;
}

void Engine_GPU::UpdateHostMirror()
{
	if (m_StaleVolt.stale)
		VoltagesToHost();
	if (m_StaleCurr.stale)
		CurrentsToHost();
}

bool Engine_GPU::IterateTS(unsigned int iterTS)
{
	// at the first call: a sub-grid is attached after Init() (Operator_CylinderMultiGrid::CreateEngine())
	if (!m_BlockingPrepared)
	{
		PrepareBlocking();
		m_BlockingPrepared = true;
	}
	TrialClock(true);
	unsigned int done = 0;
	while (done<iterTS)
	{
		// a block of k timesteps (temporal blocking), or one flat timestep
		unsigned int k = 0;
		if (m_Block>=0)
		{
			const std::vector<Engine_BlockedPlan>& plans = m_BlockCandidates[m_Block];
			k = BlockLength(iterTS-done, plans.front().k);
			if (k>=2)
				RunBlock(plans[plans.front().k - k]);
		}
		if (k<2)
		{
			k = 1;
			// an extension that needs the whole domain between the main updates (the steady-state
			// check reads the energy of all fields) gets the two main updates
			bool fused = m_Fused && !m_FlatSplit;
			for (size_t n=0; n<m_GPU_exts.size() && fused; ++n)
				fused = !m_GPU_exts[n]->NeedsFlatStep(numTS);
			if (fused)
				FusedStep();
			else
			{
				VoltageHalfStep();
				CurrentHalfStep();
			}
			NextTimestep();
		}
		done += k;
		TrialProgress(k);
	}
	TrialClock(false);
	FinishBatch();
	return true;
}

// A sub-grid follows in every half step, then the coupling of the two (the order
// of Engine_ThreadPool: the innermost coupling first).

void Engine_GPU::VoltageHalfStep()
{
	DeviceHalfStep(false, Domain(), numTS, [this]() {m_Backend->UpdateVoltages(Domain());});

	if (m_SubGrid)
	{
		static_cast<Engine_GPU*>(m_SubGrid)->VoltageHalfStep();
		SyncSubGridVoltages();
	}
}

void Engine_GPU::CurrentHalfStep()
{
	DeviceHalfStep(true, Domain(), numTS, [this]() {m_Backend->UpdateCurrents(Domain());});

	if (m_SubGrid)
	{
		static_cast<Engine_GPU*>(m_SubGrid)->CurrentHalfStep();
		SyncSubGridCurrents();
	}
}

bool Engine_GPU::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy) const
{
	return m_Backend->CalcFastEnergy(numNodes, E_energy, H_energy);
}

void Engine_GPU::NextTimestep()
{
	++numTS;
	if (m_SubGrid)
		static_cast<Engine_GPU*>(m_SubGrid)->NextTimestep();
}

void Engine_GPU::FinishBatch()
{
	// the sub-grid first fills itself from its own sub-grids (see Engine_ThreadPool::IterateTS())
	if (m_SubGrid)
	{
		static_cast<Engine_GPU*>(m_SubGrid)->FinishBatch();
		InterpolateFromSubGrid();
	}
	if (m_SharedMemory)
	{
		// the host mirror is the device memory
		VoltagesToHost();
		CurrentsToHost();
	}
	else
	{
		// read on demand, see GetVolt()
		MarkStale(m_StaleVolt);
		MarkStale(m_StaleCurr);
		// IterateTS() returns with its timesteps done, as on the CPU, e.g. for the
		// speed of the run: the kernels run asynchronously
		m_Backend->Synchronize();
	}
}
