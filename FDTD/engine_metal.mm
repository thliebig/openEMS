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

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <stdexcept>
#include <unistd.h>
#include <vector>

using std::cout;
using std::endl;

namespace
{
const char* voltageKernelSource = R"METAL(
#include <metal_stdlib>
using namespace metal;
constant bool compressedCoefficients [[function_constant(0)]];

struct GridParams
{
	uint nx;
	uint ny;
	uint nzv;
	uint start_x;
	uint num_x;
};

kernel void update_voltages(
	device float4* volt [[buffer(0)]],
	const device float4* curr [[buffer(1)]],
	const device float4* vv [[buffer(2)]],
	const device float4* vi [[buffer(3)]],
	constant GridParams& p [[buffer(4)]],
	const device uint* coeffIndex [[buffer(5), function_constant(compressedCoefficients)]],
	uint3 gid [[thread_position_in_grid]])
{
	// Packed Z is the contiguous dimension and maps to adjacent GPU lanes.
	if (gid.x >= p.nzv || gid.y >= p.ny || gid.z >= p.num_x)
		return;

	const uint z = gid.x;
	const uint y = gid.y;
	const uint x = p.start_x + gid.z;
	const uint y_stride = p.nzv * 3;
	const uint x_stride = p.ny * y_stride;
	const uint base = x * x_stride + y * y_stride + z * 3;
	const uint base_xm = x == 0 ? base : base - x_stride;
	const uint base_ym = y == 0 ? base : base - y_stride;

	const float4 cx = curr[base];
	const float4 cy = curr[base + 1];
	const float4 cz = curr[base + 2];
	const float4 hz_y = curr[base_ym + 2];
	const float4 hx_y = curr[base_ym];
	const float4 hz_x = curr[base_xm + 2];
	const float4 hy_x = curr[base_xm + 1];

	float4 hy_z;
	float4 hx_z;
	if (z == 0)
	{
		const uint end = base + (p.nzv - 1) * 3;
		const float4 hy_end = curr[end + 1];
		const float4 hx_end = curr[end];
		hy_z = float4(0.0f, hy_end.x, hy_end.y, hy_end.z);
		hx_z = float4(0.0f, hx_end.x, hx_end.y, hx_end.z);
	}
	else
	{
		hy_z = curr[base - 2];
		hx_z = curr[base - 3];
	}

	const float4 ex = volt[base];
	const float4 ey = volt[base + 1];
	const float4 ez = volt[base + 2];
	const uint c = compressedCoefficients ? coeffIndex[base / 3] * 3 : base;
	volt[base] = ex * vv[c] + vi[c] * (cz - hz_y - cy + hy_z);
	volt[base + 1] = ey * vv[c + 1] + vi[c + 1] * (cx - hx_z - cz + hz_x);
	volt[base + 2] = ez * vv[c + 2] + vi[c + 2] * (cy - hy_x - cx + hx_y);
}

kernel void update_currents(
	device float4* curr [[buffer(0)]],
	const device float4* volt [[buffer(1)]],
	const device float4* ii [[buffer(2)]],
	const device float4* iv [[buffer(3)]],
	constant GridParams& p [[buffer(4)]],
	const device uint* coeffIndex [[buffer(5), function_constant(compressedCoefficients)]],
	uint3 gid [[thread_position_in_grid]])
{
	if (gid.x >= p.nzv || gid.y >= p.ny - 1 || gid.z >= p.num_x)
		return;

	const uint z = gid.x;
	const uint y = gid.y;
	const uint x = p.start_x + gid.z;
	const uint y_stride = p.nzv * 3;
	const uint x_stride = p.ny * y_stride;
	const uint base = x * x_stride + y * y_stride + z * 3;

	const float4 ex = volt[base];
	const float4 ey = volt[base + 1];
	const float4 ez = volt[base + 2];
	const float4 ez_y = volt[base + y_stride + 2];
	const float4 ex_y = volt[base + y_stride];
	const float4 ez_x = volt[base + x_stride + 2];
	const float4 ey_x = volt[base + x_stride + 1];

	float4 ey_z;
	float4 ex_z;
	if (z + 1 < p.nzv)
	{
		ey_z = volt[base + 4];
		ex_z = volt[base + 3];
	}
	else
	{
		const uint start = base - z * 3;
		const float4 ey_start = volt[start + 1];
		const float4 ex_start = volt[start];
		ey_z = float4(ey_start.y, ey_start.z, ey_start.w, 0.0f);
		ex_z = float4(ex_start.y, ex_start.z, ex_start.w, 0.0f);
	}

	const float4 hx = curr[base];
	const float4 hy = curr[base + 1];
	const float4 hz = curr[base + 2];
	const uint c = compressedCoefficients ? coeffIndex[base / 3] * 3 : base;
	curr[base] = hx * ii[c] + iv[c] * (ez - ez_y - ey + ey_z);
	curr[base + 1] = hy * ii[c + 1] + iv[c + 1] * (ex - ex_z - ez + ez_x);
	curr[base + 2] = hz * ii[c + 2] + iv[c + 2] * (ey - ey_x - ex + ex_y);
}

// UPML arrays are scalar NIJK, while fields use the SSE z-lane layout.
// Scalar stores touch only physical cells, never SIMD padding lanes.
struct PMLParams
{
	uint sx, sy, sz;
	uint nx, ny, nz;
	uint grid_ny, grid_nzv;
};

kernel void upml_pre(
	device float* field [[buffer(0)]],
	device float* flux [[buffer(1)]],
	const device float* self [[buffer(2)]],
	const device float* oldFlux [[buffer(3)]],
	constant PMLParams& p [[buffer(4)]],
	uint gid [[thread_position_in_grid]])
{
	const uint cells = p.nx * p.ny * p.nz;
	if (gid >= 3 * cells) return;
	const uint n = gid / cells;
	const uint z = gid % p.nz + p.sz;
	const uint y = (gid / p.nz) % p.ny + p.sy;
	const uint x = (gid % cells) / (p.ny * p.nz) + p.sx;
	const uint f = (((x * p.grid_ny + y) * p.grid_nzv + z % p.grid_nzv) * 3 + n) * 4 + z / p.grid_nzv;
	const float saved = self[gid] * field[f] - oldFlux[gid] * flux[gid];
	field[f] = flux[gid];
	flux[gid] = saved;
}

kernel void upml_post(
	device float* field [[buffer(0)]],
	device float* flux [[buffer(1)]],
	const device float* newFlux [[buffer(2)]],
	constant PMLParams& p [[buffer(4)]],
	uint gid [[thread_position_in_grid]])
{
	const uint cells = p.nx * p.ny * p.nz;
	if (gid >= 3 * cells) return;
	const uint n = gid / cells;
	const uint z = gid % p.nz + p.sz;
	const uint y = (gid / p.nz) % p.ny + p.sy;
	const uint x = (gid % cells) / (p.ny * p.nz) + p.sx;
	const uint f = (((x * p.grid_ny + y) * p.grid_nzv + z % p.grid_nzv) * 3 + n) * 4 + z / p.grid_nzv;
	const float saved = flux[gid];
	flux[gid] = field[f];
	field[f] = saved + newFlux[gid] * flux[gid];
}
)METAL";

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
	id<MTLComputePipelineState> voltagePipeline;
	id<MTLComputePipelineState> currentPipeline;
	id<MTLBuffer> volt;
	id<MTLBuffer> curr;
	id<MTLBuffer> vv;
	id<MTLBuffer> vi;
	id<MTLBuffer> ii;
	id<MTLBuffer> iv;
	id<MTLBuffer> coeffIndex;
	id<MTLComputePipelineState> pmlPrePipeline;
	id<MTLComputePipelineState> pmlPostPipeline;
	id<MTLCommandBuffer> pending;

	struct PMLRegion
	{
		Engine_Ext_UPML* extension;
		PMLParams params;
		id<MTLBuffer> flux[2];
		id<MTLBuffer> self[2];
		id<MTLBuffer> oldFlux[2];
		id<MTLBuffer> newFlux[2];
	};
	std::vector<PMLRegion> pml;

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
		const char* setting = std::getenv("OPENEMS_METAL_COMPRESS");
		if (setting && setting[0] == '0')
			return;
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
		// Bound setup memory and dictionary cache footprint; require substantial
		// reuse. Highly nonuniform meshes quickly fall back to dense reads.
		const size_t maxRecords = std::min<size_t>(4096, count / 4);
		if (!maxRecords)
			return;
		std::unordered_map<Record, uint32_t, Hash> lookup;
		lookup.reserve(maxRecords);
		std::vector<uint32_t> indices(count);
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
				const uint32_t index = static_cast<uint32_t>(lookup.size());
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
				return;
		}
		id<MTLBuffer> indexBuffer = [device newBufferWithBytes:indices.data()
		    length:indices.size() * sizeof(uint32_t) options:MTLResourceStorageModeShared];
		if (!indexBuffer)
			return;
		vv = packed[0]; vi = packed[1]; ii = packed[2]; iv = packed[3];
		coeffIndex = indexBuffer;
		cout << "Metal: lossless coefficients: " << lookup.size() << " / " << count
		     << " unique packed records, " << lookup.size() * 192 + count * 4
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
	e->Init();
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
		NSString* source = [NSString stringWithUTF8String:voltageKernelSource];
		MTLCompileOptions* options = [MTLCompileOptions new];
#ifdef OPENEMS_METAL_FAST_MATH
		options.fastMathEnabled = YES;
#else
		options.fastMathEnabled = NO;
#endif
		id<MTLLibrary> library = [m_Metal->device newLibraryWithSource:source options:options error:&error];
		if (!library)
			throw MetalError("Metal: failed to compile voltage kernel", error);

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

		function = [library newFunctionWithName:@"upml_pre"];
		m_Metal->pmlPrePipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->pmlPrePipeline)
			throw MetalError("Metal: failed to create UPML pre-update pipeline", error);
		function = [library newFunctionWithName:@"upml_post"];
		m_Metal->pmlPostPipeline = [m_Metal->device newComputePipelineStateWithFunction:function error:&error];
		if (!m_Metal->pmlPostPipeline)
			throw MetalError("Metal: failed to create UPML post-update pipeline", error);
		InitUPML();

		const char* reference = std::getenv("OPENEMS_METAL_FP64_REFERENCE");
		m_Metal->referenceEnabled = reference && reference[0] != '\0' && reference[0] != '0';
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
		region.flux[0] = wrap(pml->volt_flux);
		region.flux[1] = wrap(pml->curr_flux);
		region.self[0] = wrap(op->vv);
		region.self[1] = wrap(op->ii);
		region.oldFlux[0] = wrap(op->vvfo);
		region.oldFlux[1] = wrap(op->iifo);
		region.newFlux[0] = wrap(op->vvfn);
		region.newFlux[1] = wrap(op->iifn);
		m_Metal->pml.push_back(region);
	}
	if (!m_Metal->pml.empty())
		cout << "Metal: GPU UPML conditioning: " << m_Metal->pml.size() << " regions" << endl;
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
			auto region = std::find_if(m_Metal->pml.begin(), m_Metal->pml.end(),
				[extension](const MetalState::PMLRegion& r) { return r.extension == extension; });
			if (region == m_Metal->pml.end())
			{
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
			id<MTLComputePipelineState> pipeline = pre ? m_Metal->pmlPrePipeline : m_Metal->pmlPostPipeline;
			id<MTLComputeCommandEncoder> encoder = [m_Metal->Commands() computeCommandEncoder];
			[encoder setComputePipelineState:pipeline];
			[encoder setBuffer:voltage ? m_Metal->volt : m_Metal->curr offset:0 atIndex:0];
			[encoder setBuffer:region->flux[f] offset:0 atIndex:1];
			[encoder setBuffer:pre ? region->self[f] : region->newFlux[f] offset:0 atIndex:2];
			if (pre)
				[encoder setBuffer:region->oldFlux[f] offset:0 atIndex:3];
			[encoder setBytes:&region->params length:sizeof(PMLParams) atIndex:4];
			const PMLParams& p = region->params;
			[encoder dispatchThreads:MTLSizeMake((NSUInteger)p.nx * p.ny * p.nz * 3, 1, 1)
				threadsPerThreadgroup:MTLSizeMake(pipeline.threadExecutionWidth, 1, 1)];
			[encoder endEncoding];
		}
		if (!pre)
			FinishMetalCommands(); // Apply hooks, sources, probes and dumps run on CPU.
	}
}

void Engine_Metal::DoPreVoltageUpdates() { RunUPMLExtensions(true, true); }
void Engine_Metal::DoPostVoltageUpdates() { RunUPMLExtensions(true, false); }
void Engine_Metal::DoPreCurrentUpdates() { RunUPMLExtensions(false, true); }
void Engine_Metal::DoPostCurrentUpdates() { RunUPMLExtensions(false, false); }

void Engine_Metal::Reset()
{
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
		if (m_Metal->referenceEnabled || m_Metal->pml.empty())
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
		if (m_Metal->referenceEnabled || m_Metal->pml.empty())
			FinishMetalCommands();
	}
	if (m_Metal->referenceEnabled)
	{
		m_Metal->CompareAndSnapshot(m_Metal->refCurr, m_Metal->lastCurr,
		                            f4_curr_ptr->data(), f4_curr_ptr->size()*4, false);
	}
}
