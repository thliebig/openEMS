/*
* Copyright (C) 2026 openEMS contributors
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*/

#include "engine_metal.h"
#include "extensions/engine_ext_upml.h"
#include "extensions/operator_ext_upml.h"
#include "extensions/engine_ext_excitation.h"
#include "extensions/operator_ext_excitation.h"
#include "extensions/engine_ext_lorentzmaterial.h"
#include "excitation.h"
#include "metal_library.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdexcept>
#include <unistd.h>
#include <utility>
#include <vector>

using std::cout;
using std::endl;

namespace
{

struct PMLParams
{
	uint32_t sx, sy, sz;
	uint32_t nx, ny, nz;
	uint32_t grid_ny, grid_nzv;
};

struct GridParams
{
	uint32_t nx;
	uint32_t ny;
	uint32_t nzv;
	uint32_t start_x;
	uint32_t num_x;
};

struct ExcitationSource
{
	uint32_t fieldIndex;
	float amplitude;
	uint32_t delay;
};

struct ExcitationParams
{
	uint32_t count;
	uint32_t timestep;
	uint32_t signalLength;
	uint32_t period;
};

static const unsigned int DIAMOND_DEPTH = 4;

struct DiamondStep
{
	int32_t voltageRange[4];
	int32_t currentRange[4];
	uint32_t voltageSourceOffset, voltageSourceCount;
	uint32_t currentSourceOffset, currentSourceCount;
};

struct DiamondTile { DiamondStep steps[DIAMOND_DEPTH]; };
struct DiamondParams { uint32_t nx, ny, nzv, timestep, depth; };

struct DiamondSource
{
	uint32_t fieldIndex;
	float amplitude;
	uint32_t delay, signalOffset, signalLength, period;
};

using DiamondRange = std::pair<int32_t, int32_t>;
using DiamondBlock = std::vector<DiamondRange>;
using DiamondAxis = std::array<std::vector<DiamondBlock>, 2>;

// Mountain/valley construction adapted for Metal from the experimental
// project-diamond-rework1 CPU tiler. At every half-step, both phases together
// partition the axis while each phase remains internally independent.
static DiamondAxis MakeDiamondAxis(uint32_t width, uint32_t blockWidth, uint32_t halfSteps)
{
	const int shortest = blockWidth;
	const int longest = blockWidth + halfSteps - 1;
	int blocks = width / (shortest + longest) * 2;
	int remainder = width % (shortest + longest);
	for (int n = 0; remainder > 0; ++n)
	{
		++blocks;
		remainder -= n % 2 == 0 ? shortest : longest;
	}
	std::vector<DiamondBlock> all(blocks + 1, DiamondBlock(halfSteps));
	int last = -1;
	for (size_t n = 0; n < all.size(); ++n)
	{
		const int span = n % 2 == 0 ? shortest : longest;
		all[n][halfSteps - 1] = {last + 1, last + span};
		last += span;
	}
	for (int t = halfSteps - 2; t >= 0; --t)
		for (size_t n = 0; n < all.size(); ++n)
		{
			const DiamondRange next = all[n][t + 1];
			DiamondRange range;
			if (n % 2 == 0)
				range = t % 2 ? DiamondRange(next.first - 1, next.second)
				              : DiamondRange(next.first, next.second + 1);
			else
				range = t % 2 ? DiamondRange(next.first, next.second - 1)
				              : DiamondRange(next.first + 1, next.second);
			range.first = std::max<int32_t>(range.first, 0);
			all[n][t] = range;
		}
	DiamondAxis phases;
	for (size_t n = 0; n < all.size(); ++n)
	{
		for (DiamondRange& range : all[n])
		{
			if (range.first >= (int32_t)width)
				range = {-1, -1};
			else
				range.second = std::min<int32_t>(range.second, width - 1);
		}
		phases[n % 2].push_back(all[n]);
	}
	return phases;
}

static_assert(sizeof(DiamondStep) == 48, "Metal DiamondStep ABI");
static_assert(sizeof(DiamondTile) == 192, "Metal DiamondTile ABI");
static_assert(sizeof(DiamondParams) == 20, "Metal DiamondParams ABI");
static_assert(sizeof(DiamondSource) == 24, "Metal DiamondSource ABI");

static std::runtime_error MetalError(const char* what, NSError* error)
{
	std::string message(what);
	if (error)
	{
		message += ": ";
		message += [[error localizedDescription] UTF8String];
	}
	return std::runtime_error(message);
}
}

struct Engine_Metal::MetalState
{
	id<MTLDevice> device;
	id<MTLCommandQueue> queue;
	id<MTLLibrary> library;
	id<MTLComputePipelineState> voltagePipeline;
	id<MTLComputePipelineState> currentPipeline;
	id<MTLComputePipelineState> diamondPipeline;
	id<MTLComputePipelineState> excitationPipeline;
	id<MTLComputePipelineState> adeAdvancePipeline;
	id<MTLComputePipelineState> adeApplyPipeline;
	id<MTLBuffer> volt;
	id<MTLBuffer> curr;
	id<MTLBuffer> vv;
	id<MTLBuffer> vi;
	id<MTLBuffer> ii;
	id<MTLBuffer> iv;
	id<MTLBuffer> coeffIndex;
	id<MTLComputePipelineState> pmlPrePipeline[2];
	id<MTLComputePipelineState> pmlPostPipeline[2];
	size_t reusedPMLBytes = 0;
	id<MTLCommandBuffer> pending;
	bool diamondRequested = true;
	bool legacyRequested = false;
	id<MTLBuffer> diamondTiles[DIAMOND_DEPTH + 1][4];
	id<MTLBuffer> diamondSourceIndices[DIAMOND_DEPTH + 1][4];
	id<MTLBuffer> diamondSourceTable;
	uint32_t diamondTileCount[DIAMOND_DEPTH + 1][4] = {};
	bool diamondUpdate = false;
	bool diamondHasSources = false;
	id<MTLBuffer> diamondSignal;

	struct ExcitationRegion
	{
		Engine_Ext_Excitation* extension;
		id<MTLBuffer> sources[2];
		id<MTLBuffer> signals[2];
		uint32_t counts[2];
		uint32_t signalLength;
		uint32_t period;
	};
	std::vector<ExcitationRegion> excitations;

	struct PMLRegion
	{
		Engine_Ext_UPML* extension;
		PMLParams params;
		id<MTLBuffer> indices;
		id<MTLBuffer> coeffIndex[2];
		id<MTLBuffer> flux[2];
		id<MTLBuffer> self[2];
		id<MTLBuffer> oldFlux[2];
		id<MTLBuffer> newFlux[2];
		bool coefficientsReleased[2] = {false, false};
	};
	std::vector<PMLRegion> pml;

	// Plain volt-ADE (conducting-sheet) regions offloaded to the GPU. One entry
	// per active packed-field edge, with two poles packed into coeff/state.
	struct ADERegion
	{
		Engine_Ext_LorentzMaterial* extension;
		uint32_t count;
		id<MTLBuffer> indices;
		id<MTLBuffer> coeff;
		id<MTLBuffer> state;
	};
	std::vector<ADERegion> ade;

	void CompressPML(PMLRegion& region, unsigned int field)
	{
		// Lossless dictionaries substantially reduce resident PML memory.
		using Record = std::array<uint32_t, 3>;
		struct Hash { size_t operator()(const Record& r) const {
			return ((size_t)r[0]*16777619U ^ r[1])*16777619U ^ r[2];
		}};
		const size_t count = (size_t)region.params.nx * region.params.ny * region.params.nz * 3;
		const uint32_t* arrays[] = {static_cast<const uint32_t*>(region.self[field].contents),
			static_cast<const uint32_t*>(region.oldFlux[field].contents),
			static_cast<const uint32_t*>(region.newFlux[field].contents)};
		std::unordered_map<Record, uint16_t, Hash> lookup;
		std::vector<uint16_t> indices(count);
		std::vector<uint32_t> records[3];
		const size_t limit = 65536;
		for (size_t i=0; i<count; ++i)
		{
			Record record = {{arrays[0][i], arrays[1][i], arrays[2][i]}};
			auto found = lookup.find(record);
			if (found == lookup.end())
			{
				if (lookup.size() >= limit) {
					cout << "Metal: UPML coefficient dictionary limit reached; using dense UPML coefficients" << endl;
					return;
				}
				uint16_t id = static_cast<uint16_t>(lookup.size());
				found = lookup.emplace(record, id).first;
				for (size_t a=0; a<3; ++a) records[a].push_back(record[a]);
			}
			indices[i] = found->second;
		}
		id<MTLBuffer> packed[3];
		for (size_t a=0; a<3; ++a) {
			packed[a] = [device newBufferWithBytes:records[a].data() length:records[a].size()*4 options:MTLResourceStorageModeShared];
			if (!packed[a]) {
				cout << "Metal: UPML dictionary buffer allocation failed; using dense UPML coefficients" << endl;
				return;
			}
		}
		id<MTLBuffer> index = [device newBufferWithBytes:indices.data() length:indices.size()*2 options:MTLResourceStorageModeShared];
		if (!index) {
			cout << "Metal: UPML index buffer allocation failed; using dense UPML coefficients" << endl;
			return;
		}
		region.self[field]=packed[0]; region.oldFlux[field]=packed[1]; region.newFlux[field]=packed[2];
		region.coeffIndex[field]=index;
		region.coefficientsReleased[field]=true;
		ArrayLib::ArrayNIJK<FDTD_FLOAT>* dense[] = {
			field ? &region.extension->m_Op_UPML->ii : &region.extension->m_Op_UPML->vv,
			field ? &region.extension->m_Op_UPML->iifo : &region.extension->m_Op_UPML->vvfo,
			field ? &region.extension->m_Op_UPML->iifn : &region.extension->m_Op_UPML->vvfn};
		for (auto* array : dense)
		{
			reordered.erase(std::remove_if(reordered.begin(), reordered.end(),
				[array](const ReorderedArray& saved) { return saved.array == array; }), reordered.end());
			array->Reset();
		}
		const size_t compactBytes=count*2+lookup.size()*12;
		cout << "Metal: lossless UPML coefficients: " << lookup.size() << "/" << count
		     << " records, " << compactBytes << " GPU bytes; released " << count*12
		     << " dense bytes" << endl;
	}

	// Coefficient owners are the operator and may outlive/recreate the engine.
	// Restore their scalar order on Reset; engine-only flux arrays die with it.
	struct ReorderedArray
	{
		ArrayLib::ArrayNIJK<FDTD_FLOAT>* array;
		PMLParams params;
		id<MTLBuffer> indices;
	};
	std::vector<ReorderedArray> reordered;
	// Allocate before mutating operator-owned arrays; teardown must not allocate.
	std::vector<float> reorderScratch;

	static size_t ScalarIndex(const PMLParams& p, uint32_t fieldIndex)
	{
		uint32_t lane = fieldIndex % 4;
		fieldIndex /= 4;
		uint32_t component = fieldIndex % 3;
		fieldIndex /= 3;
		uint32_t z = fieldIndex % p.grid_nzv + lane * p.grid_nzv;
		fieldIndex /= p.grid_nzv;
		uint32_t y = fieldIndex % p.grid_ny, x = fieldIndex / p.grid_ny;
		const size_t cells = (size_t)p.nx * p.ny * p.nz;
		return component * cells + ((size_t)(x-p.sx)*p.ny + y-p.sy)*p.nz + z-p.sz;
	}

	void RestorePML()
	{
		for (const auto& saved : reordered)
		{
			const PMLParams& p = saved.params;
			const uint32_t* indices = static_cast<const uint32_t*>(saved.indices.contents);
			float* scalar = reorderScratch.data();
			for (size_t i = 0; i < saved.array->size(); ++i)
				scalar[ScalarIndex(p, indices[i])] = saved.array->data()[i];
			std::memcpy(saved.array->data(), scalar, saved.array->bytes());
		}
		reordered.clear();

		// Compact dictionaries replaced and released the operator-owned dense
		// arrays during the run. Recreate scalar NIJK order so the same operator
		// remains valid if a caller constructs another CPU or Metal engine.
		for (auto& region : pml)
			for (unsigned int field=0; field<2; ++field)
			{
				if (!region.coefficientsReleased[field]) continue;
				const PMLParams& p=region.params;
				const size_t count=(size_t)p.nx*p.ny*p.nz*3;
				const uint32_t* fieldIndices=static_cast<const uint32_t*>(region.indices.contents);
				const uint16_t* dictionaryIndices=static_cast<const uint16_t*>(region.coeffIndex[field].contents);
				ArrayLib::ArrayNIJK<FDTD_FLOAT>* arrays[] = {
					field ? &region.extension->m_Op_UPML->ii : &region.extension->m_Op_UPML->vv,
					field ? &region.extension->m_Op_UPML->iifo : &region.extension->m_Op_UPML->vvfo,
					field ? &region.extension->m_Op_UPML->iifn : &region.extension->m_Op_UPML->vvfn};
				id<MTLBuffer> buffers[] = {region.self[field],region.oldFlux[field],region.newFlux[field]};
				for (unsigned int a=0; a<3; ++a)
				{
					arrays[a]->Init("restored_upml", {p.nx,p.ny,p.nz});
					const FDTD_FLOAT* dictionary=static_cast<const FDTD_FLOAT*>(buffers[a].contents);
					for (size_t i=0; i<count; ++i)
						arrays[a]->data()[ScalarIndex(p,fieldIndices[i])]=dictionary[dictionaryIndices[i]];
				}
				region.coefficientsReleased[field]=false;
			}
		reorderScratch.clear();
		reorderScratch.shrink_to_fit();
	}

	id<MTLCommandBuffer> Commands()
	{
		if (!pending)
			pending = [queue commandBuffer];
		if (!pending)
			throw std::runtime_error("Metal: failed to create command buffer");
		return pending;
	}

	// Coefficients are immutable during stepping. Keep the operator's dense
	// arrays for CPU access and the FP64 reference; only GPU reads use this copy.
	void CompressCoefficients(size_t count)
	{
		using Record = std::array<uint32_t, 48>; // 4 arrays * 3 components * 4 lanes
		struct Hash
		{
			size_t operator()(const Record& record) const
			{
				size_t hash = 14695981039346656037ULL;
				for (uint32_t word : record)
					hash = (hash ^ word) * 1099511628211ULL;
				return hash;
			}
		};
		// Bound setup memory and dictionary cache footprint; require at least
		// fourfold reuse. The packed index is a uint16, so 65536 records
		// (indices 0..65535) is the format limit.
		size_t maxRecords = std::min<size_t>(65536, count / 4);
		if (!maxRecords)
		{
			cout << "Metal: too few positions to compress; using dense operator coefficients" << endl;
			return;
		}
		std::unordered_map<Record, uint16_t, Hash> lookup;
		lookup.reserve(maxRecords);
		std::vector<uint16_t> indices(count);
		std::vector<uint32_t> dictionaries[4];
		const void* dense[] = {vv.contents, vi.contents, ii.contents, iv.contents};
		Record previous{};
		for (size_t pos = 0; pos < count; ++pos)
		{
			Record record;
			for (size_t a = 0; a < 4; ++a)
				std::memcpy(record.data() + a * 12,
				    static_cast<const char*>(dense[a]) + pos * 48, 48);
			// Uniform runs need neither hashing nor a dictionary search.
			if (pos && record == previous)
			{
				indices[pos] = indices[pos - 1];
				continue;
			}
			previous = record;
			auto found = lookup.find(record);
			if (found == lookup.end())
			{
				if (lookup.size() == maxRecords)
				{
					cout << "Metal: coefficient dictionary limit reached; using dense coefficients" << endl;
					return;
				}
				const uint16_t index = static_cast<uint16_t>(lookup.size());
				found = lookup.emplace(record, index).first;
				for (size_t a = 0; a < 4; ++a)
					dictionaries[a].insert(dictionaries[a].end(),
					    record.begin() + a * 12, record.begin() + (a + 1) * 12);
			}
			indices[pos] = found->second;
		}
		id<MTLBuffer> packed[4];
		for (size_t a = 0; a < 4; ++a)
		{
			packed[a] = [device newBufferWithBytes:dictionaries[a].data()
			    length:dictionaries[a].size() * sizeof(uint32_t) options:MTLResourceStorageModeShared];
			if (!packed[a])
			{
				cout << "Metal: coefficient dictionary buffer allocation failed; using dense coefficients" << endl;
				return;
			}
		}
		id<MTLBuffer> indexBuffer = [device newBufferWithBytes:indices.data()
		    length:indices.size() * sizeof(uint16_t) options:MTLResourceStorageModeShared];
		if (!indexBuffer)
		{
			cout << "Metal: coefficient index buffer allocation failed; using dense coefficients" << endl;
			return;
		}
		vv = packed[0]; vi = packed[1]; ii = packed[2]; iv = packed[3];
		coeffIndex = indexBuffer;
		cout << "Metal: lossless coefficients: " << lookup.size() << " / " << count
		     << " unique packed records, " << lookup.size() * 192 + count * 2
		     << " GPU bytes (dense " << count * 192 << ")" << endl;
	}

	bool hasUPML = false;
	bool referenceEnabled = false;
	std::vector<double> refVolt;
	std::vector<double> refCurr;
	std::vector<float> lastVolt;
	std::vector<float> lastCurr;
	double voltageSquaredError = 0;
	double voltageSquaredReference = 0;
	double currentSquaredError = 0;
	double currentSquaredReference = 0;
	double voltageMaxAbs = 0;
	double currentMaxAbs = 0;
	uint64_t voltageSamples = 0;
	uint64_t currentSamples = 0;

	void Reconcile(std::vector<double>& ref, std::vector<float>& last,
	               const f4vector* actual, size_t count)
	{
		const float* values = reinterpret_cast<const float*>(actual);
		for (size_t n=0; n<count; ++n)
		{
			// Flux swaps are not additive field corrections. With UPML, check
			// each stencil from the actual conditioned input instead of evolving
			// an unstable reference that has no matching FP64 flux state.
			if (hasUPML)
				ref[n] = values[n];
			else
				ref[n] += (double)values[n] - (double)last[n];
			last[n] = values[n];
		}
	}

	void CompareAndSnapshot(const std::vector<double>& ref, std::vector<float>& last,
	                        const f4vector* actual, size_t count, bool voltage)
	{
		const float* values = reinterpret_cast<const float*>(actual);
		double squaredError = 0;
		double squaredReference = 0;
		double maxAbs = 0;
		for (size_t n=0; n<count; ++n)
		{
			double error = std::abs((double)values[n] - ref[n]);
			maxAbs = std::max(maxAbs, error);
			squaredError += error * error;
			squaredReference += ref[n] * ref[n];
			last[n] = values[n];
		}
		if (voltage)
		{
			voltageSquaredError += squaredError;
			voltageSquaredReference += squaredReference;
			voltageMaxAbs = std::max(voltageMaxAbs, maxAbs);
			voltageSamples += count;
		}
		else
		{
			currentSquaredError += squaredError;
			currentSquaredReference += squaredReference;
			currentMaxAbs = std::max(currentMaxAbs, maxAbs);
			currentSamples += count;
		}
	}
};

Engine_Metal* Engine_Metal::New(const Operator_sse* op)
{
	cout << "Create FDTD engine (Metal field updates)" << endl;
	Engine_Metal* e = new Engine_Metal(op);
	try { e->Init(); }
	catch (...) { delete e; throw; }
	return e;
}

Engine_Metal::Engine_Metal(const Operator_sse* op) : Engine_sse(op)
{
	m_Metal = NULL;
}

Engine_Metal::~Engine_Metal()
{
	Reset();
}

void Engine_Metal::Init()
{
	Engine_sse::Init();

	@autoreleasepool
	{
		m_Metal = new MetalState();
		m_Metal->device = MTLCreateSystemDefaultDevice();
		if (!m_Metal->device)
			throw std::runtime_error("Metal: no GPU device available");

		m_Metal->queue = [m_Metal->device newCommandQueue];
		if (!m_Metal->queue)
			throw std::runtime_error("Metal: failed to create command queue");

		NSError* error = nil;
		m_Metal->library = OpenEMSMetalLibrary(m_Metal->device, &error);
		if (!m_Metal->library)
			throw MetalError("Metal: failed to load precompiled kernel library", error);
		id<MTLLibrary> library = m_Metal->library;

		const NSUInteger pageSize = (NSUInteger)getpagesize();
		auto paddedLength = [pageSize](NSUInteger bytes) {
			return (bytes + pageSize - 1) / pageSize * pageSize;
		};
		auto noFree = ^(void*, NSUInteger) {};
		const NSUInteger fieldBytes = f4_volt_ptr->bytes();
		const NSUInteger coeffBytes = Op->f4_vv_ptr->bytes();
		m_Metal->volt = [m_Metal->device newBufferWithBytesNoCopy:f4_volt_ptr->data()
			length:paddedLength(fieldBytes) options:MTLResourceStorageModeShared deallocator:noFree];
		m_Metal->curr = [m_Metal->device newBufferWithBytesNoCopy:f4_curr_ptr->data()
			length:paddedLength(fieldBytes) options:MTLResourceStorageModeShared deallocator:noFree];
		m_Metal->vv = [m_Metal->device newBufferWithBytesNoCopy:Op->f4_vv_ptr->data()
			length:paddedLength(coeffBytes) options:MTLResourceStorageModeShared deallocator:noFree];
		m_Metal->vi = [m_Metal->device newBufferWithBytesNoCopy:Op->f4_vi_ptr->data()
			length:paddedLength(coeffBytes) options:MTLResourceStorageModeShared deallocator:noFree];
		m_Metal->ii = [m_Metal->device newBufferWithBytesNoCopy:Op->f4_ii_ptr->data()
			length:paddedLength(coeffBytes) options:MTLResourceStorageModeShared deallocator:noFree];
		m_Metal->iv = [m_Metal->device newBufferWithBytesNoCopy:Op->f4_iv_ptr->data()
			length:paddedLength(coeffBytes) options:MTLResourceStorageModeShared deallocator:noFree];
		if (!m_Metal->volt || !m_Metal->curr || !m_Metal->vv || !m_Metal->vi ||
			!m_Metal->ii || !m_Metal->iv)
			throw std::runtime_error("Metal: failed to wrap shared buffers");

		m_Metal->CompressCoefficients(f4_volt_ptr->size() / 3);
		MTLFunctionConstantValues* constants = [MTLFunctionConstantValues new];
		bool compressed = m_Metal->coeffIndex != nil;
		[constants setConstantValue:&compressed type:MTLDataTypeBool atIndex:0];
		id<MTLFunction> function = [library newFunctionWithName:@"update_voltages"
		    constantValues:constants error:&error];
		if (!function)
			throw MetalError("Metal: update_voltages kernel not found", error);
		m_Metal->voltagePipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->voltagePipeline)
			throw MetalError("Metal: failed to create voltage pipeline", error);
		function = [library newFunctionWithName:@"update_currents" constantValues:constants error:&error];
		if (!function)
			throw MetalError("Metal: update_currents kernel not found", error);
		m_Metal->currentPipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->currentPipeline)
			throw MetalError("Metal: failed to create current pipeline", error);
		function = [library newFunctionWithName:@"apply_excitation"];
		if (!function)
			throw std::runtime_error("Metal: apply_excitation kernel not found");
		m_Metal->excitationPipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->excitationPipeline)
			throw MetalError("Metal: failed to create excitation pipeline", error);
		function = [library newFunctionWithName:@"ade_advance"];
		if (!function)
			throw std::runtime_error("Metal: ade_advance kernel not found");
		m_Metal->adeAdvancePipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->adeAdvancePipeline)
			throw MetalError("Metal: failed to create ADE advance pipeline", error);
		function = [library newFunctionWithName:@"ade_apply"];
		if (!function)
			throw std::runtime_error("Metal: ade_apply kernel not found");
		m_Metal->adeApplyPipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->adeApplyPipeline)
			throw MetalError("Metal: failed to create ADE apply pipeline", error);

		for (unsigned int compressed=0; compressed<2; ++compressed)
		{
			MTLFunctionConstantValues* values = [MTLFunctionConstantValues new];
			bool enabled = compressed != 0;
			[values setConstantValue:&enabled type:MTLDataTypeBool atIndex:1];
			function = [library newFunctionWithName:@"upml_indexed_pre" constantValues:values error:&error];
			if (!function) throw MetalError("Metal: UPML pre kernel specialization failed", error);
			m_Metal->pmlPrePipeline[compressed] = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
			if (!m_Metal->pmlPrePipeline[compressed]) throw MetalError("Metal: UPML pre pipeline failed", error);
			function = [library newFunctionWithName:@"upml_indexed_post" constantValues:values error:&error];
			if (!function) throw MetalError("Metal: UPML post kernel specialization failed", error);
			m_Metal->pmlPostPipeline[compressed] = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
			if (!m_Metal->pmlPostPipeline[compressed]) throw MetalError("Metal: UPML post pipeline failed", error);
		}
		InitUPML();
		InitExcitations();

		const char* reference = std::getenv("OPENEMS_METAL_FP64_REFERENCE");
		m_Metal->referenceEnabled = reference && reference[0] != '\0' && reference[0] != '0';
		const char* wavefront = std::getenv("OPENEMS_METAL_FUSED_PIPELINE");
		if (wavefront && wavefront[0] == '0')
		{
			m_Metal->diamondRequested = false;
			m_Metal->legacyRequested = true;
			std::cerr << "Metal: legacy two-dispatch pipeline selected by OPENEMS_METAL_FUSED_PIPELINE=0" << std::endl;
		}
		if (m_Metal->referenceEnabled)
		{
			m_Metal->diamondRequested = false;
			m_Metal->legacyRequested = true;
		}
		InitADE();
		if (!m_Metal->ade.empty())
			m_Metal->diamondRequested = false;
		InitDiamondUpdate();
		cout << "Metal: in-place diamond E/H pipeline: " << (m_Metal->diamondUpdate ? "enabled" : "disabled") << endl;
		if (m_Metal->referenceEnabled)
		{
			size_t scalarCount = f4_volt_ptr->size() * 4;
			m_Metal->refVolt.assign(scalarCount, 0.0);
			m_Metal->refCurr.assign(scalarCount, 0.0);
			m_Metal->lastVolt.assign(scalarCount, 0.0f);
			m_Metal->lastCurr.assign(scalarCount, 0.0f);
			cout << "Metal: enabled diagnostic FP64 update reference" << endl;
			if (m_Metal->hasUPML)
				cout << "Metal: UPML FP64 reference checks local stencils, not flux evolution" << endl;
		}
	}
}

void Engine_Metal::InitUPML()
{
	for (Engine_Extension* extension : m_Eng_exts)
		if (dynamic_cast<Engine_Ext_UPML*>(extension))
			m_Metal->hasUPML = true;
	const char* setting = std::getenv("OPENEMS_METAL_PML");
	if (setting && setting[0] == '0')
	{
		m_Metal->diamondRequested = false;
		m_Metal->legacyRequested = true;
		cout << "Metal: CPU UPML conditioning selected" << endl;
		return;
	}

	const NSUInteger pageSize = (NSUInteger)getpagesize();
	auto wrap = [&](const ArrayLib::ArrayNIJK<FDTD_FLOAT>& array) {
		const NSUInteger bytes = (array.bytes() + pageSize - 1) / pageSize * pageSize;
		id<MTLBuffer> buffer = [m_Metal->device newBufferWithBytesNoCopy:array.data()
			length:bytes options:MTLResourceStorageModeShared deallocator:^(void*, NSUInteger) {}];
		if (!buffer)
			throw std::runtime_error("Metal: failed to wrap UPML buffer");
		return buffer;
	};
	for (Engine_Extension* extension : m_Eng_exts)
	{
		Engine_Ext_UPML* pml = dynamic_cast<Engine_Ext_UPML*>(extension);
		if (!pml)
			continue;
		Operator_Ext_UPML* op = pml->m_Op_UPML;
		// Opposing slabs can leave an empty interior for another face.
		// Keep its no-op CPU hook rather than creating zero-length buffers.
		if (!op->m_numLines[0] || !op->m_numLines[1] || !op->m_numLines[2])
			continue;
		MetalState::PMLRegion region;
		region.extension = pml;
		region.params = {op->m_StartPos[0], op->m_StartPos[1], op->m_StartPos[2],
			op->m_numLines[0], op->m_numLines[1], op->m_numLines[2], numLines[1], numVectors};
		{
			const PMLParams& p = region.params;
			// The grid index range was validated in Operator_Metal::SetupCSXGrid.
			const size_t cells = (size_t)p.nx * p.ny * p.nz;
			// Precompute the permutation once. Reuse CPU storage in indexed order;
			// CPU UPML hooks never execute on these buffers during Metal stepping.
			const size_t componentCount=cells*3;
			region.indices = [m_Metal->device newBufferWithLength:componentCount * sizeof(uint32_t)
				options:MTLResourceStorageModeShared];
			if (!region.indices)
				throw std::runtime_error("Metal: failed to allocate UPML indices");
			uint32_t* indices = static_cast<uint32_t*>(region.indices.contents);
			uint32_t* const indicesBegin=indices;
			for (uint32_t x = 0; x < p.nx; ++x)
				for (uint32_t y = 0; y < p.ny; ++y)
					for (uint32_t zv = 0; zv < p.grid_nzv; ++zv)
						for (uint32_t n = 0; n < 3; ++n)
							for (uint32_t lane = 0; lane < 4; ++lane)
							{
								const uint32_t z = zv + lane * p.grid_nzv;
								if (z < p.sz || z >= p.sz + p.nz) continue;
								const size_t f = ((((size_t)(x + p.sx) * p.grid_ny + y + p.sy)
									* p.grid_nzv + zv) * 3 + n) * 4 + lane;
								*indices++ = static_cast<uint32_t>(f);
							}
			if ((size_t)(indices-indicesBegin) != componentCount)
				throw std::runtime_error("Metal: incomplete UPML index mapping");
			if (m_Metal->reorderScratch.size() < componentCount)
				m_Metal->reorderScratch.resize(componentCount);
			auto pack = [&](ArrayLib::ArrayNIJK<FDTD_FLOAT>& array, bool restore) {
				// Reuse the operator-owned CPU storage in indexed order; CPU UPML
				// hooks never execute on it during Metal stepping.
				id<MTLBuffer> buffer = wrap(array);
				float* out = m_Metal->reorderScratch.data();
				// Register before mutation so exception cleanup can restore the operator.
				if (restore) m_Metal->reordered.push_back({&array, p, region.indices});
				const uint32_t* fieldIndices=static_cast<const uint32_t*>(region.indices.contents);
				for (size_t i = 0; i < componentCount; ++i)
					out[i] = array.data()[MetalState::ScalarIndex(p,fieldIndices[i])];
				std::memcpy(array.data(), out, array.bytes());
				m_Metal->reusedPMLBytes += array.bytes();
				return buffer;
			};
			region.flux[0] = pack(pml->volt_flux, false);
			region.flux[1] = pack(pml->curr_flux, false);
			region.self[0] = pack(op->vv, true);
			region.self[1] = pack(op->ii, true);
			region.oldFlux[0] = pack(op->vvfo, true);
			region.oldFlux[1] = pack(op->iifo, true);
			region.newFlux[0] = pack(op->vvfn, true);
			region.newFlux[1] = pack(op->iifn, true);
		}
		// Register before compressing: CompressPML may release the operator-owned
		// dense arrays, and RestorePML can only rebuild them from pml.
		m_Metal->pml.push_back(region);
		m_Metal->CompressPML(m_Metal->pml.back(), 0);
		m_Metal->CompressPML(m_Metal->pml.back(), 1);
	}
	// When every reordered operator coefficient was replaced by a compact
	// dictionary, no scalar restoration needs the setup scratch. Flux storage is
	// engine-owned and intentionally remains indexed until engine destruction.
	if (m_Metal->reordered.empty())
	{
		m_Metal->reorderScratch.clear();
		m_Metal->reorderScratch.shrink_to_fit();
	}
	if (!m_Metal->pml.empty())
	{
		cout << "Metal: GPU UPML conditioning: " << m_Metal->pml.size() << " regions" << endl;
		cout << "Metal: UPML duplicate bytes avoided: " << m_Metal->reusedPMLBytes << endl;
		cout << "Metal: UPML reorder scratch bytes retained: " << m_Metal->reorderScratch.size()*sizeof(float) << endl;
	}
}

void Engine_Metal::InitExcitations()
{
	for (Engine_Extension* extension : m_Eng_exts)
	{
		Engine_Ext_Excitation* excitation = dynamic_cast<Engine_Ext_Excitation*>(extension);
		if (!excitation)
		{
			if (!dynamic_cast<Engine_Ext_UPML*>(extension) && m_Metal->diamondRequested)
			{
				m_Metal->diamondRequested = false;
				std::cerr << "Metal: extension '" << extension->GetExtensionName()
				          << "' has not migrated to the diamond wavefront" << std::endl;
			}
			continue;
		}
		Operator_Ext_Excitation* op = excitation->m_Op_Exc;
		if (!op || !op->m_Exc || op->m_Exc->GetLength() == 0)
		{
			if (m_Metal->diamondRequested)
			{
				m_Metal->diamondRequested = false;
				std::cerr << "Metal: excitation has no signal and cannot use the diamond wavefront" << std::endl;
			}
			continue;
		}
		MetalState::ExcitationRegion region;
		region.extension = excitation;
		region.counts[0] = op->Volt_Count;
		region.counts[1] = op->Curr_Count;
		region.signalLength = op->m_Exc->GetLength();
		double signalPeriod = op->m_Exc->GetSignalPeriod();
		region.period = signalPeriod > 0 ? static_cast<uint32_t>(signalPeriod / op->m_Exc->GetTimestep()) : 0;
		for (unsigned int field = 0; field < 2; ++field)
		{
			const uint32_t count = region.counts[field];
			if (!count) continue;
			std::vector<ExcitationSource> sources(count);
			for (uint32_t n = 0; n < count; ++n)
			{
				const uint32_t x = field ? op->Curr_index[0][n] : op->Volt_index[0][n];
				const uint32_t y = field ? op->Curr_index[1][n] : op->Volt_index[1][n];
				const uint32_t z = field ? op->Curr_index[2][n] : op->Volt_index[2][n];
				const uint32_t dir = field ? op->Curr_dir[n] : op->Volt_dir[n];
				const size_t index = ((((size_t)x * numLines[1] + y) * numVectors + z % numVectors) * 3 + dir) * 4 + z / numVectors;
				sources[n] = {static_cast<uint32_t>(index),
					field ? op->Curr_amp[n] : op->Volt_amp[n],
					field ? op->Curr_delay[n] : op->Volt_delay[n]};
			}
			region.sources[field] = [m_Metal->device newBufferWithBytes:sources.data()
				length:sources.size() * sizeof(ExcitationSource) options:MTLResourceStorageModeShared];
			FDTD_FLOAT* signal = field ? op->m_Exc->GetCurrentSignal() : op->m_Exc->GetVoltageSignal();
			region.signals[field] = [m_Metal->device newBufferWithBytes:signal
				length:region.signalLength * sizeof(FDTD_FLOAT) options:MTLResourceStorageModeShared];
			if (!region.sources[field] || !region.signals[field])
				throw std::runtime_error("Metal: failed to allocate excitation buffers");
		}
		m_Metal->excitations.push_back(region);
	}
}

void Engine_Metal::InitDiamondUpdate()
{
	if (!m_Metal->diamondRequested || !m_Metal->pml.empty() || !m_Metal->ade.empty())
	{
		if (m_Metal->legacyRequested)
			return;
		const std::string reason = !m_Metal->pml.empty() ? "UPML" :
			!m_Metal->ade.empty() ? "ADE" : "a CPU extension hook";
		throw std::runtime_error("Metal: " + reason +
			" has not migrated to the diamond E/H kernel; refusing to start a legacy simulation. "
			"Use OPENEMS_METAL_FUSED_PIPELINE=0 only for an explicit diagnostic run.");
	}
	struct SourceTemplate { uint32_t sourceIndex, x, y; bool voltage; };
	std::vector<SourceTemplate> sourceTemplates;
	std::vector<DiamondSource> sourceTable;
	std::vector<float> signal;
	for (const auto& region : m_Metal->excitations)
		for (unsigned int field = 0; field < 2; ++field)
		{
			if (!region.counts[field])
				continue;
			if (signal.size() + region.signalLength > std::numeric_limits<uint32_t>::max())
				throw std::runtime_error("Metal: diamond signal buffer exceeds uint32");
			const uint32_t signalOffset = static_cast<uint32_t>(signal.size());
			const float* samples = static_cast<const float*>(region.signals[field].contents);
			signal.insert(signal.end(), samples, samples + region.signalLength);
			const ExcitationSource* sources = static_cast<const ExcitationSource*>(region.sources[field].contents);
			for (uint32_t n = 0; n < region.counts[field]; ++n)
			{
				if (sourceTable.size() == std::numeric_limits<uint32_t>::max())
					throw std::runtime_error("Metal: diamond source table exceeds uint32");
				uint32_t q = sources[n].fieldIndex / 4 / 3 / numVectors;
				const uint32_t y = q % numLines[1];
				const uint32_t x = q / numLines[1];
				sourceTemplates.push_back({static_cast<uint32_t>(sourceTable.size()),
					x, y, field == 0});
				sourceTable.push_back({sources[n].fieldIndex, sources[n].amplitude,
					sources[n].delay, signalOffset, region.signalLength, region.period});
			}
		}
	NSError* error = nil;
	MTLFunctionConstantValues* constants = [MTLFunctionConstantValues new];
	bool compressed = m_Metal->coeffIndex != nil;
	bool hasSources = !sourceTemplates.empty();
	[constants setConstantValue:&compressed type:MTLDataTypeBool atIndex:0];
	[constants setConstantValue:&hasSources type:MTLDataTypeBool atIndex:2];
	id<MTLFunction> function = [m_Metal->library newFunctionWithName:@"update_diamond"
		constantValues:constants error:&error];
	if (!function)
		throw MetalError("Metal: update_diamond kernel not found", error);
	m_Metal->diamondPipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
	if (!m_Metal->diamondPipeline)
		throw MetalError("Metal: failed to create diamond pipeline", error);
	m_Metal->diamondHasSources = hasSources;
	if (hasSources)
	{
		m_Metal->diamondSourceTable = [m_Metal->device newBufferWithBytes:sourceTable.data()
			length:sourceTable.size() * sizeof(DiamondSource) options:MTLResourceStorageModeShared];
		m_Metal->diamondSignal = [m_Metal->device newBufferWithBytes:signal.data()
			length:signal.size() * sizeof(float) options:MTLResourceStorageModeShared];
		if (!m_Metal->diamondSourceTable || !m_Metal->diamondSignal)
			throw std::runtime_error("Metal: failed to allocate diamond source buffers");
	}

	std::array<std::vector<std::vector<uint32_t>>, 2> sourcesByXY = {
		std::vector<std::vector<uint32_t>>((size_t)numLines[0] * numLines[1]),
		std::vector<std::vector<uint32_t>>((size_t)numLines[0] * numLines[1])};
	for (const SourceTemplate& source : sourceTemplates)
		sourcesByXY[source.voltage ? 0 : 1][(size_t)source.x * numLines[1] + source.y]
			.push_back(source.sourceIndex);

	const uint32_t blockWidth = 2;
	uint64_t auxiliaryBytes = signal.size() * sizeof(float) +
		sourceTable.size() * sizeof(DiamondSource);
	for (uint32_t depth = 1; depth <= DIAMOND_DEPTH; ++depth)
	{
		const DiamondAxis xAxis = MakeDiamondAxis(numLines[0], blockWidth, depth * 2);
		const DiamondAxis yAxis = MakeDiamondAxis(numLines[1], blockWidth, depth * 2);
		for (uint32_t phase = 0; phase < 4; ++phase)
		{
			const uint32_t phaseX = phase / 2, phaseY = phase % 2;
			std::vector<DiamondTile> tiles;
			for (const DiamondBlock& xb : xAxis[phaseX])
				for (const DiamondBlock& yb : yAxis[phaseY])
				{
					DiamondTile tile{};
					bool any = false;
					for (uint32_t t = 0; t < depth; ++t)
					{
						const DiamondRange ex = xb[2 * t], hx = xb[2 * t + 1];
						const DiamondRange ey = yb[2 * t], hy = yb[2 * t + 1];
						DiamondStep& step = tile.steps[t];
						for (int n = 0; n < 4; ++n)
							step.voltageRange[n] = step.currentRange[n] = -1;
						if (ex.first >= 0 && ey.first >= 0)
						{
							step.voltageRange[0] = ex.first; step.voltageRange[1] = ex.second;
							step.voltageRange[2] = ey.first; step.voltageRange[3] = ey.second;
							any = true;
						}
						if (hx.first >= 0 && hy.first >= 0)
						{
							step.currentRange[0] = hx.first; step.currentRange[1] = hx.second;
							step.currentRange[2] = hy.first; step.currentRange[3] = hy.second;
							any = true;
						}
					}
					if (any)
						tiles.push_back(tile);
				}
			std::vector<uint32_t> sourceIndices;
			auto appendSource = [&sourceIndices](uint32_t source) {
				if (sourceIndices.size() == std::numeric_limits<uint32_t>::max())
					throw std::runtime_error("Metal: diamond source schedule exceeds uint32");
				sourceIndices.push_back(source);
			};
			auto appendRange = [&](const int32_t range[4], unsigned int field) {
				if (range[0] < 0)
					return;
				for (int32_t x = range[0]; x <= range[1]; ++x)
					for (int32_t y = range[2]; y <= range[3]; ++y)
						for (uint32_t source : sourcesByXY[field][(size_t)x * numLines[1] + y])
							appendSource(source);
			};
			for (DiamondTile& tile : tiles)
				for (uint32_t t = 0; t < depth; ++t)
				{
					DiamondStep& step = tile.steps[t];
					step.voltageSourceOffset = static_cast<uint32_t>(sourceIndices.size());
					appendRange(step.voltageRange, 0);
					step.voltageSourceCount = static_cast<uint32_t>(sourceIndices.size()) - step.voltageSourceOffset;
					step.currentSourceOffset = static_cast<uint32_t>(sourceIndices.size());
					appendRange(step.currentRange, 1);
					step.currentSourceCount = static_cast<uint32_t>(sourceIndices.size()) - step.currentSourceOffset;
				}
			if (tiles.empty())
				continue;
			auxiliaryBytes += tiles.size() * sizeof(DiamondTile) +
				sourceIndices.size() * sizeof(uint32_t);
			m_Metal->diamondTileCount[depth][phase] = static_cast<uint32_t>(tiles.size());
			m_Metal->diamondTiles[depth][phase] = [m_Metal->device newBufferWithBytes:tiles.data()
				length:tiles.size() * sizeof(DiamondTile) options:MTLResourceStorageModeShared];
			if (hasSources)
				m_Metal->diamondSourceIndices[depth][phase] = sourceIndices.empty() ?
					[m_Metal->device newBufferWithLength:sizeof(uint32_t) options:MTLResourceStorageModeShared] :
					[m_Metal->device newBufferWithBytes:sourceIndices.data()
					 length:sourceIndices.size() * sizeof(uint32_t) options:MTLResourceStorageModeShared];
			if (!m_Metal->diamondTiles[depth][phase] ||
				(hasSources && !m_Metal->diamondSourceIndices[depth][phase]))
				throw std::runtime_error("Metal: failed to allocate diamond schedule");
		}
	}
	m_Metal->diamondUpdate = true;
	cout << "Metal: in-place diamond update: " << DIAMOND_DEPTH
	     << " timesteps/block, width " << blockWidth
	     << ", " << auxiliaryBytes << " auxiliary bytes" << endl;
}

void Engine_Metal::ApplyMetalExcitations(bool voltage)
{
	@autoreleasepool
	{
		const unsigned int field = voltage ? 0 : 1;
		for (const auto& region : m_Metal->excitations)
		{
			if (!region.counts[field]) continue;
			ExcitationParams params = {region.counts[field], numTS, region.signalLength,
				region.period ? region.period : numTS + 1};
			id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
			[encoder setComputePipelineState:m_Metal->excitationPipeline];
			[encoder setBuffer:voltage ? m_Metal->volt : m_Metal->curr offset:0 atIndex:0];
			[encoder setBuffer:region.sources[field] offset:0 atIndex:1];
			[encoder setBuffer:region.signals[field] offset:0 atIndex:2];
			[encoder setBytes:&params length:sizeof(params) atIndex:3];
			[encoder dispatchThreads:MTLSizeMake(1, 1, 1) threadsPerThreadgroup:MTLSizeMake(1, 1, 1)];
			[encoder endEncoding];
		}
	}
}

void Engine_Metal::FinishMetalCommands()
{
	if (!m_Metal || !m_Metal->pending)
		return;
	id<MTLCommandBuffer> commands = m_Metal->pending;
	[commands commit];
	[commands waitUntilCompleted];
	m_Metal->pending = nil;
	if (commands.status == MTLCommandBufferStatusError)
		throw MetalError("Metal: field/UPML update failed", commands.error);
}

void Engine_Metal::RunUPMLExtensions(bool voltage, bool pre)
{
	@autoreleasepool
	{
		// Keep the CPU extension order exactly: pre in reverse, post forward.
		// Separate encoders order even overlapping regions; CPU hooks see only
		// completed GPU writes. Typically pre/field/post share one submission.
		for (size_t i = 0; i < m_Eng_exts.size(); ++i)
		{
			Engine_Extension* extension = m_Eng_exts[pre ? m_Eng_exts.size() - 1 - i : i];
			// Offloaded ADE extensions advance before the voltage update, never
			// drain the GPU, and never run the CPU hook. Their apply pass runs from
			// Apply2Voltages after all post updates.
			if (HasADEOffload(extension))
			{
				if (voltage && pre)
					AdvanceADEOffload(extension);
				continue;
			}
			auto region = std::find_if(m_Metal->pml.begin(), m_Metal->pml.end(),
				[extension](const MetalState::PMLRegion& r) { return r.extension == extension; });
			if (region == m_Metal->pml.end())
			{
				if (m_Metal->diamondRequested && dynamic_cast<Engine_Ext_Excitation*>(extension))
					continue;
				FinishMetalCommands();
				if (voltage)
				{
					if (pre) extension->DoPreVoltageUpdates();
					else extension->DoPostVoltageUpdates();
				}
				else
				{
					if (pre) extension->DoPreCurrentUpdates();
					else extension->DoPostCurrentUpdates();
				}
				continue;
			}
			const unsigned int f = voltage ? 0 : 1;
			unsigned int compressed = region->coeffIndex[f] != nil;
			id<MTLComputePipelineState> pipeline = pre ? m_Metal->pmlPrePipeline[compressed] : m_Metal->pmlPostPipeline[compressed];
			id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
			[encoder setComputePipelineState:pipeline];
			if (compressed) [encoder setBuffer:region->coeffIndex[f] offset:0 atIndex:6];
			[encoder setBuffer:voltage ? m_Metal->volt : m_Metal->curr offset:0 atIndex:0];
			[encoder setBuffer:region->flux[f] offset:0 atIndex:1];
			[encoder setBuffer:pre ? region->self[f] : region->newFlux[f] offset:0 atIndex:2];
			if (pre)
				[encoder setBuffer:region->oldFlux[f] offset:0 atIndex:3];
			[encoder setBuffer:region->indices offset:0 atIndex:5];
			const PMLParams& p = region->params;
			[encoder dispatchThreads:MTLSizeMake((NSUInteger)p.nx * p.ny * p.nz * 3, 1, 1)
				threadsPerThreadgroup:MTLSizeMake(pipeline.threadExecutionWidth, 1, 1)];
			[encoder endEncoding];
		}
		if (!pre && !m_Metal->diamondRequested)
			FinishMetalCommands(); // CPU Apply hooks need completed fields.
	}
}

void Engine_Metal::DoPreVoltageUpdates() { RunUPMLExtensions(true, true); }
void Engine_Metal::DoPostVoltageUpdates() { RunUPMLExtensions(true, false); }
void Engine_Metal::DoPreCurrentUpdates() { RunUPMLExtensions(false, true); }
void Engine_Metal::DoPostCurrentUpdates() { RunUPMLExtensions(false, false); }

bool Engine_Metal::HasADEOffload(const Engine_Extension* extension) const
{
	for (const auto& region : m_Metal->ade)
		if (region.extension == extension)
			return true;
	return false;
}

void Engine_Metal::AdvanceADEOffload(Engine_Extension* extension)
{
	for (auto& region : m_Metal->ade)
	{
		if (region.extension != extension)
			continue;
		id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
		[encoder setComputePipelineState:m_Metal->adeAdvancePipeline];
		[encoder setBuffer:m_Metal->volt offset:0 atIndex:0];
		[encoder setBuffer:region.state offset:0 atIndex:1];
		[encoder setBuffer:region.coeff offset:0 atIndex:2];
		[encoder setBuffer:region.indices offset:0 atIndex:3];
		[encoder setBytes:&region.count length:sizeof(region.count) atIndex:4];
		[encoder dispatchThreads:MTLSizeMake((NSUInteger)region.count, 1, 1)
			threadsPerThreadgroup:MTLSizeMake(m_Metal->adeAdvancePipeline.threadExecutionWidth, 1, 1)];
		[encoder endEncoding];
		return;
	}
}

void Engine_Metal::ApplyADEOffload(Engine_Extension* extension)
{
	for (auto& region : m_Metal->ade)
	{
		if (region.extension != extension)
			continue;
		id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
		[encoder setComputePipelineState:m_Metal->adeApplyPipeline];
		[encoder setBuffer:m_Metal->volt offset:0 atIndex:0];
		[encoder setBuffer:region.state offset:0 atIndex:1];
		[encoder setBuffer:region.indices offset:0 atIndex:2];
		[encoder setBytes:&region.count length:sizeof(region.count) atIndex:3];
		[encoder dispatchThreads:MTLSizeMake((NSUInteger)region.count, 1, 1)
			threadsPerThreadgroup:MTLSizeMake(m_Metal->adeApplyPipeline.threadExecutionWidth, 1, 1)];
		[encoder endEncoding];
		return;
	}
}

void Engine_Metal::Apply2Voltages()
{
	@autoreleasepool
	{
		// CPU apply hooks first (the base class order), then the GPU ADE apply.
		// ADE edges live outside every PML region, so their relative order to the
		// other extensions cannot alias the same field entries.
		for (Engine_Extension* extension : m_Eng_exts)
			if (!HasADEOffload(extension))
				extension->Apply2Voltages();
		for (Engine_Extension* extension : m_Eng_exts)
			ApplyADEOffload(extension);
	}
}

void Engine_Metal::InitADE()
{
	// The FP64 reference reproduces every extension on the CPU, and the diamond
	// pipeline never calls Apply2Voltages, so both keep the ADE on the CPU.
	if (m_Metal->referenceEnabled || m_Metal->diamondRequested)
	{
		for (Engine_Extension* extension : m_Eng_exts)
		{
			Engine_Ext_LorentzMaterial* lor = dynamic_cast<Engine_Ext_LorentzMaterial*>(extension);
			if (lor && lor->MetalADEOffloadSupported())
			{
				std::cerr << "Metal: conducting-sheet ADE stays on the CPU ("
				          << (m_Metal->referenceEnabled ? "FP64 reference mode" : "diamond pipeline")
				          << ")" << std::endl;
				break;
			}
		}
		return;
	}

	const uint32_t ny = numLines[1];
	const uint32_t nzv = numVectors;
	size_t total = 0;
	for (Engine_Extension* extension : m_Eng_exts)
	{
		Engine_Ext_LorentzMaterial* lor = dynamic_cast<Engine_Ext_LorentzMaterial*>(extension);
		if (!lor)
			continue;
		if (!lor->MetalADEOffloadSupported())
		{
			// Lorentz flux states and ADE-current schemes need extra GPU state.
			std::cerr << "Metal: dispersive/ADE extension '" << extension->GetExtensionName()
			          << "' stays on the CPU (only the conducting-sheet volt-ADE is offloaded)" << std::endl;
			continue;
		}

		// In the conducting-sheet model both ADE poles share one field position,
		// so merge by packed field index. One thread then owns every pole of one
		// edge, which keeps the apply subtraction race-free and bit-identical to
		// the CPU sequence (pole 0 then pole 1).
		std::unordered_map<uint32_t, uint32_t> first;
		std::vector<uint32_t> indices;
		std::vector<float> coeff;
		const int order = lor->MetalADEOrder();
		for (int o = 0; o < order; ++o)
		{
			if (!lor->MetalADEVoltOn(o))
				continue;
			const unsigned int count = lor->MetalADECount(o);
			const unsigned int* px = lor->MetalADEPos(o, 0);
			const unsigned int* py = lor->MetalADEPos(o, 1);
			const unsigned int* pz = lor->MetalADEPos(o, 2);
			for (unsigned int i = 0; i < count; ++i)
			{
				const uint32_t x = px[i], y = py[i], z = pz[i];
				const uint32_t slot = z % nzv;
				const uint32_t lane = z / nzv;
				for (int n = 0; n < 3; ++n)
				{
					const float vi = lor->MetalADEVoltInt(o, n)[i];
					const float ve = lor->MetalADEVoltExt(o, n)[i];
					if (vi == 0.0f && ve == 0.0f)
						continue;
					const uint32_t f = (3 * ((x * ny + y) * nzv + slot) + n) * 4 + lane;
					auto found = first.find(f);
					uint32_t entry;
					if (found == first.end())
					{
						entry = static_cast<uint32_t>(indices.size());
						first.emplace(f, entry);
						indices.push_back(f);
						// Identity poles so an absent order stays a no-op.
						coeff.insert(coeff.end(), {1.0f, 0.0f, 1.0f, 0.0f});
					}
					else
						entry = found->second;
					coeff[entry * 4 + 2 * o] = vi;
					coeff[entry * 4 + 2 * o + 1] = ve;
				}
			}
		}
		if (indices.empty())
			continue;

		const NSUInteger stateBytes = indices.size() * 2 * sizeof(float);
		MetalState::ADERegion region;
		region.extension = lor;
		region.count = static_cast<uint32_t>(indices.size());
		region.indices = [m_Metal->device newBufferWithBytes:indices.data()
			length:indices.size() * sizeof(uint32_t) options:MTLResourceStorageModeShared];
		region.coeff = [m_Metal->device newBufferWithBytes:coeff.data()
			length:coeff.size() * sizeof(float) options:MTLResourceStorageModeShared];
		region.state = [m_Metal->device newBufferWithLength:stateBytes options:MTLResourceStorageModeShared];
		if (!region.indices || !region.coeff || !region.state)
			throw std::runtime_error("Metal: failed to allocate ADE buffers");
		std::memset(region.state.contents, 0, stateBytes);
		total += indices.size();
		m_Metal->ade.push_back(region);
	}
	if (!m_Metal->ade.empty())
		cout << "Metal: ADE offload: " << m_Metal->ade.size() << " region(s), "
		     << total << " active edges" << endl;
}

void Engine_Metal::UpdateDiamond(unsigned int depth)
{
	DiamondParams params = {numLines[0], numLines[1], numVectors, numTS, depth};
	// The kernel maps whole packed-Z slot groups to threads, so the threadgroup
	// must be an exact multiple of the slot count it can cover.
	const NSUInteger maxThreads = std::min<NSUInteger>({
		m_Metal->diamondPipeline.maxTotalThreadsPerThreadgroup,
		m_Metal->device.maxThreadsPerThreadgroup.width});
	const NSUInteger slots = std::min<NSUInteger>(numVectors, maxThreads);
	const NSUInteger groups = std::max<NSUInteger>(1, 256 / slots);
	const NSUInteger threads = slots * groups;
	for (unsigned int phase = 0; phase < 4; ++phase)
	{
		const uint32_t count = m_Metal->diamondTileCount[depth][phase];
		if (!count)
			continue;
		id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
		[encoder setComputePipelineState:m_Metal->diamondPipeline];
		[encoder setBuffer:m_Metal->volt offset:0 atIndex:0];
		[encoder setBuffer:m_Metal->curr offset:0 atIndex:1];
		[encoder setBuffer:m_Metal->vv offset:0 atIndex:2];
		[encoder setBuffer:m_Metal->vi offset:0 atIndex:3];
		[encoder setBuffer:m_Metal->ii offset:0 atIndex:4];
		[encoder setBuffer:m_Metal->iv offset:0 atIndex:5];
		[encoder setBytes:&params length:sizeof(params) atIndex:6];
		[encoder setBuffer:m_Metal->diamondTiles[depth][phase] offset:0 atIndex:7];
		if (m_Metal->coeffIndex)
			[encoder setBuffer:m_Metal->coeffIndex offset:0 atIndex:8];
		if (m_Metal->diamondHasSources)
		{
			[encoder setBuffer:m_Metal->diamondSourceTable offset:0 atIndex:9];
			[encoder setBuffer:m_Metal->diamondSourceIndices[depth][phase] offset:0 atIndex:10];
			[encoder setBuffer:m_Metal->diamondSignal offset:0 atIndex:11];
		}
		[encoder dispatchThreadgroups:MTLSizeMake(count, 1, 1)
			threadsPerThreadgroup:MTLSizeMake(threads, 1, 1)];
		[encoder endEncoding];
	}
}

bool Engine_Metal::IterateTS(unsigned int iterTS)
{
	if (m_Metal->diamondUpdate)
	{
		unsigned int remaining = iterTS;
		while (remaining)
		{
			const unsigned int depth = std::min<unsigned int>(DIAMOND_DEPTH, remaining);
			UpdateDiamond(depth);
			numTS += depth;
			remaining -= depth;
		}
		FinishMetalCommands();
		return true;
	}
	if (!m_Metal->diamondRequested)
		return Engine::IterateTS(iterTS);

	for (unsigned int iter = 0; iter < iterTS; ++iter)
	{
		DoPreVoltageUpdates();
		UpdateVoltages(0, numLines[0]);
		DoPostVoltageUpdates();
		if (m_Metal->referenceEnabled)
		{
			FinishMetalCommands();
			Engine::Apply2Voltages();
		}
		else
			ApplyMetalExcitations(true);

		DoPreCurrentUpdates();
		UpdateCurrents(0, numLines[0] - 1);
		DoPostCurrentUpdates();
		if (m_Metal->referenceEnabled)
		{
			FinishMetalCommands();
			Engine::Apply2Current();
		}
		else
			ApplyMetalExcitations(false);

		// Submit one ordered GPU pipeline per timestep. This is also the
		// synchronization point for CPU probes, dumps, and convergence checks.
		FinishMetalCommands();
		++numTS;
	}
	return true;
}

void Engine_Metal::Reset()
{
	// RestorePML rebuilds operator-owned coefficient arrays and therefore
	// allocates. Reset runs from the destructor, so a failure must not escape.
	if (m_Metal)
	{
		try { FinishMetalCommands(); }
		catch (const std::exception& e) { std::cerr << "Metal: failed to finish commands on reset: " << e.what() << std::endl; }
		catch (...) { std::cerr << "Metal: failed to finish commands on reset" << std::endl; }
		try { m_Metal->RestorePML(); }
		catch (const std::exception& e) { std::cerr << "Metal: failed to restore UPML coefficients on reset: " << e.what() << std::endl; }
		catch (...) { std::cerr << "Metal: failed to restore UPML coefficients on reset" << std::endl; }
	}
	if (m_Metal && m_Metal->referenceEnabled)
	{
		double voltL2 = m_Metal->voltageSquaredReference > 0 ?
			std::sqrt(m_Metal->voltageSquaredError / m_Metal->voltageSquaredReference) : 0;
		double currL2 = m_Metal->currentSquaredReference > 0 ?
			std::sqrt(m_Metal->currentSquaredError / m_Metal->currentSquaredReference) : 0;
		cout << std::setprecision(8)
		     << "Metal FP64 update reference: E max abs=" << m_Metal->voltageMaxAbs
		     << ", relative L2=" << voltL2
		     << "; H max abs=" << m_Metal->currentMaxAbs
		     << ", relative L2=" << currL2 << endl;
	}
	delete m_Metal;
	m_Metal = NULL;
	Engine_sse::Reset();
}

void Engine_Metal::UpdateVoltages(unsigned int startX, unsigned int numX)
{
	if (numX == 0)
		return;

	if (m_Metal->referenceEnabled)
	{
		FinishMetalCommands(); // The FP64 diagnostic needs the conditioned fields.
		const size_t count = f4_volt_ptr->size() * 4;
		m_Metal->Reconcile(m_Metal->refVolt, m_Metal->lastVolt, f4_volt_ptr->data(), count);
		m_Metal->Reconcile(m_Metal->refCurr, m_Metal->lastCurr, f4_curr_ptr->data(), count);
		const float* vv = reinterpret_cast<const float*>(Op->f4_vv_ptr->data());
		const float* vi = reinterpret_cast<const float*>(Op->f4_vi_ptr->data());
		const uint32_t ny = numLines[1], nzv = numVectors;
		auto idx = [ny,nzv](uint32_t n,uint32_t x,uint32_t y,uint32_t z,uint32_t l) {
			return (size_t)((((x*ny+y)*nzv+z)*3+n)*4+l);
		};
		for (uint32_t x=startX; x<startX+numX; ++x)
			for (uint32_t y=0; y<ny; ++y)
				for (uint32_t z=0; z<nzv; ++z)
					for (uint32_t l=0; l<4; ++l)
					{
						size_t v0=idx(0,x,y,z,l), v1=idx(1,x,y,z,l), v2=idx(2,x,y,z,l);
						size_t hzY=idx(2,x,y?y-1:y,z,l), hxY=idx(0,x,y?y-1:y,z,l);
						size_t hzX=idx(2,x?x-1:x,y,z,l), hyX=idx(1,x?x-1:x,y,z,l);
						double hyZ, hxZ;
						if (z) { hyZ=m_Metal->refCurr[idx(1,x,y,z-1,l)]; hxZ=m_Metal->refCurr[idx(0,x,y,z-1,l)]; }
						else if (l) { hyZ=m_Metal->refCurr[idx(1,x,y,nzv-1,l-1)]; hxZ=m_Metal->refCurr[idx(0,x,y,nzv-1,l-1)]; }
						else { hyZ=0; hxZ=0; }
						double cx=m_Metal->refCurr[v0], cy=m_Metal->refCurr[v1], cz=m_Metal->refCurr[v2];
						m_Metal->refVolt[v0]=m_Metal->refVolt[v0]*vv[v0]+vi[v0]*(cz-m_Metal->refCurr[hzY]-cy+hyZ);
						m_Metal->refVolt[v1]=m_Metal->refVolt[v1]*vv[v1]+vi[v1]*(cx-hxZ-cz+m_Metal->refCurr[hzX]);
						m_Metal->refVolt[v2]=m_Metal->refVolt[v2]*vv[v2]+vi[v2]*(cy-m_Metal->refCurr[hyX]-cx+m_Metal->refCurr[hxY]);
					}
	}

	@autoreleasepool
	{
		id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
		[encoder setComputePipelineState:m_Metal->voltagePipeline];
		[encoder setBuffer:m_Metal->volt offset:0 atIndex:0];
		[encoder setBuffer:m_Metal->curr offset:0 atIndex:1];
		[encoder setBuffer:m_Metal->vv offset:0 atIndex:2];
		[encoder setBuffer:m_Metal->vi offset:0 atIndex:3];
		if (m_Metal->coeffIndex)
			[encoder setBuffer:m_Metal->coeffIndex offset:0 atIndex:5];

		GridParams params = {numLines[0], numLines[1], numVectors, startX, numX};
		[encoder setBytes:&params length:sizeof(params) atIndex:4];

		NSUInteger width = m_Metal->voltagePipeline.threadExecutionWidth;
		MTLSize threadsPerGroup = MTLSizeMake(width, 1, 1);
		MTLSize grid = MTLSizeMake(numVectors, numLines[1], numX);
		[encoder dispatchThreads:grid threadsPerThreadgroup:threadsPerGroup];
		[encoder endEncoding];
		if (m_Metal->referenceEnabled || (!m_Metal->diamondRequested && m_Metal->pml.empty()))
			FinishMetalCommands();
	}
	if (m_Metal->referenceEnabled)
	{
		m_Metal->CompareAndSnapshot(m_Metal->refVolt, m_Metal->lastVolt,
		                            f4_volt_ptr->data(), f4_volt_ptr->size()*4, true);
	}
}

void Engine_Metal::UpdateCurrents(unsigned int startX, unsigned int numX)
{
	if (numX == 0)
		return;

	if (m_Metal->referenceEnabled)
	{
		FinishMetalCommands();
		// Reproduce voltage post-conditioning and current pre-conditioning.
		m_Metal->Reconcile(m_Metal->refVolt, m_Metal->lastVolt,
		                   f4_volt_ptr->data(), f4_volt_ptr->size()*4);
		m_Metal->Reconcile(m_Metal->refCurr, m_Metal->lastCurr,
		                   f4_curr_ptr->data(), f4_curr_ptr->size()*4);
		const float* ii = reinterpret_cast<const float*>(Op->f4_ii_ptr->data());
		const float* iv = reinterpret_cast<const float*>(Op->f4_iv_ptr->data());
		const uint32_t ny = numLines[1], nzv = numVectors;
		auto idx = [ny,nzv](uint32_t n,uint32_t x,uint32_t y,uint32_t z,uint32_t l) {
			return (size_t)((((x*ny+y)*nzv+z)*3+n)*4+l);
		};
		for (uint32_t x=startX; x<startX+numX; ++x)
			for (uint32_t y=0; y<ny-1; ++y)
				for (uint32_t z=0; z<nzv; ++z)
					for (uint32_t l=0; l<4; ++l)
					{
						size_t i0=idx(0,x,y,z,l), i1=idx(1,x,y,z,l), i2=idx(2,x,y,z,l);
						size_t ezY=idx(2,x,y+1,z,l), exY=idx(0,x,y+1,z,l);
						size_t ezX=idx(2,x+1,y,z,l), eyX=idx(1,x+1,y,z,l);
						double eyZ, exZ;
						if (z+1<nzv) { eyZ=m_Metal->refVolt[idx(1,x,y,z+1,l)]; exZ=m_Metal->refVolt[idx(0,x,y,z+1,l)]; }
						else if (l<3) { eyZ=m_Metal->refVolt[idx(1,x,y,0,l+1)]; exZ=m_Metal->refVolt[idx(0,x,y,0,l+1)]; }
						else { eyZ=0; exZ=0; }
						double ex=m_Metal->refVolt[i0], ey=m_Metal->refVolt[i1], ez=m_Metal->refVolt[i2];
						m_Metal->refCurr[i0]=m_Metal->refCurr[i0]*ii[i0]+iv[i0]*(ez-m_Metal->refVolt[ezY]-ey+eyZ);
						m_Metal->refCurr[i1]=m_Metal->refCurr[i1]*ii[i1]+iv[i1]*(ex-exZ-ez+m_Metal->refVolt[ezX]);
						m_Metal->refCurr[i2]=m_Metal->refCurr[i2]*ii[i2]+iv[i2]*(ey-m_Metal->refVolt[eyX]-ex+m_Metal->refVolt[exY]);
					}
	}

	@autoreleasepool
	{
		id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
		[encoder setComputePipelineState:m_Metal->currentPipeline];
		[encoder setBuffer:m_Metal->curr offset:0 atIndex:0];
		[encoder setBuffer:m_Metal->volt offset:0 atIndex:1];
		[encoder setBuffer:m_Metal->ii offset:0 atIndex:2];
		[encoder setBuffer:m_Metal->iv offset:0 atIndex:3];
		if (m_Metal->coeffIndex)
			[encoder setBuffer:m_Metal->coeffIndex offset:0 atIndex:5];

		GridParams params = {numLines[0], numLines[1], numVectors, startX, numX};
		[encoder setBytes:&params length:sizeof(params) atIndex:4];

		NSUInteger width = m_Metal->currentPipeline.threadExecutionWidth;
		MTLSize threadsPerGroup = MTLSizeMake(width, 1, 1);
		MTLSize grid = MTLSizeMake(numVectors, numLines[1] - 1, numX);
		[encoder dispatchThreads:grid threadsPerThreadgroup:threadsPerGroup];
		[encoder endEncoding];
		if (m_Metal->referenceEnabled || (!m_Metal->diamondRequested && m_Metal->pml.empty()))
			FinishMetalCommands();
	}
	if (m_Metal->referenceEnabled)
	{
		m_Metal->CompareAndSnapshot(m_Metal->refCurr, m_Metal->lastCurr,
		                            f4_curr_ptr->data(), f4_curr_ptr->size()*4, false);
	}
}
