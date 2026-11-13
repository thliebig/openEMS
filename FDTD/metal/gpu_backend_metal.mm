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
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "metal_internal.h"
#include "metal_cm_source.h"
#include "FDTD/operator.h"
#include "FDTD/engine_plan.h"
#include "FDTD/gpu_cell_model.h"
#include "FDTD/kernels/ext_ops.h"

/***************************** Impl *****************************/

id<MTLComputePipelineState> Metal_Context::Pipeline(const char* function)
{
	std::map<std::string, id<MTLComputePipelineState>>::iterator it = pipelines.find(function);
	if (it!=pipelines.end())
		return it->second;

	@autoreleasepool
	{
		MTLCompileOptions* options = [MTLCompileOptions new];
		// keep IEEE float semantics, the engines must agree with the CPU
#if __MAC_OS_X_VERSION_MAX_ALLOWED >= 150000
		if (@available(macOS 15.0, *))
			options.mathMode = MTLMathModeSafe;
		else
#endif
		{
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
			options.fastMathEnabled = NO;
#pragma clang diagnostic pop
		}

		// all kernels in one library, compiled once
		NSError* error = nil;
		if (!library)
		{
			library = [device newLibraryWithSource:[NSString stringWithUTF8String:METAL_CM_SOURCE] options:options error:&error];
			if (!library)
				throw std::runtime_error(std::string("GPU_Backend_Metal: kernel compilation failed: ") + [[error localizedDescription] UTF8String]);
		}
		id<MTLFunction> fn = [library newFunctionWithName:[NSString stringWithUTF8String:function]];
		if (!fn)
			throw std::runtime_error(std::string("GPU_Backend_Metal: kernel not found: ") + function);
		id<MTLComputePipelineState> pso = [device newComputePipelineStateWithFunction:fn error:&error];
		if (!pso)
			throw std::runtime_error(std::string("GPU_Backend_Metal: pipeline creation failed: ") + [[error localizedDescription] UTF8String]);
		pipelines[function] = pso;
		return pso;
	}
}

id<MTLComputeCommandEncoder> Metal_Context::Encoder()
{
	if (enc)
		return enc;
	// command buffer and encoder are autoreleased: without a pool they would pile up for a run
	@autoreleasepool
	{
		cmd = [queue commandBuffer];
		// concurrent dispatch with explicit barriers, see GPU_Backend_Metal::Impl::Dispatch()
		enc = [cmd computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
		// the buffers the kernels reach through pointers in their arguments, in one call
		if (!resources.empty())
		{
			std::vector<__unsafe_unretained id<MTLResource>> list(resources.begin(), resources.end());
			[enc useResources:list.data() count:list.size() usage:(MTLResourceUsageRead | MTLResourceUsageWrite)];
		}
	}
	return enc;
}

id<MTLCommandBuffer> Metal_Context::Commit()
{
	if (!cmd)
		return nil;
	id<MTLCommandBuffer> committed = cmd;
	@autoreleasepool
	{
		[enc endEncoding];
		[cmd commit];
		enc = nil;
		cmd = nil;
	}
	return committed;
}

void Metal_Context::Flush()
{
	if (!cmd)
		return;
	@autoreleasepool
	{
		[enc endEncoding];
		[cmd commit];
		[cmd waitUntilCompleted];
		// done with this batch either way: a failed one must not be ended or committed again
		const bool failed = (cmd.status==MTLCommandBufferStatusError);
		const std::string what = failed ? [[cmd.error localizedDescription] UTF8String] : "";
		enc = nil;
		cmd = nil;
		if (failed)
			throw std::runtime_error("GPU_Backend_Metal: command buffer failed: " + what);
	}
}

id<MTLBuffer> GPU_Backend_Metal::Impl::NewBuffer(size_t bytes, const void* data)
{
	// shared storage: the host and the kernels work on the same pages, so the fields never
	// have to be transferred. Metal rejects a zero length buffer, empty arrays get the minimum.
	id<MTLBuffer> buf = [device newBufferWithLength:std::max(bytes, (size_t)4) options:MTLResourceStorageModeShared];
	if (!buf)
		throw std::runtime_error("GPU_Backend_Metal: buffer allocation failed");
	std::memset([buf contents], 0, [buf length]);
	if (data && bytes)
		std::memcpy([buf contents], data, bytes);
	return buf;
}

uint64_t GPU_Backend_Metal::Impl::NewResource(size_t bytes, const void* data)
{
	id<MTLBuffer> buf = NewBuffer(bytes, data);
	ctx->resources.push_back(buf);
	if (ctx->enc)   // an encoder already open declares it too
		[ctx->enc useResource:buf usage:(MTLResourceUsageRead | MTLResourceUsageWrite)];
	uint64_t address = 0;
	if (@available(macOS 13.0, *))
		address = buf.gpuAddress;
	ctx->addresses[address] = buf;
	return address;
}

void GPU_Backend_Metal::Impl::Dispatch(id<MTLComputePipelineState> pso, size_t ni, size_t nj, size_t nk)
{
	if (ni==0 || nj==0 || nk==0)
		return;
	id<MTLComputeCommandEncoder> enc = Encoder();
	[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
	// wide in the first dimension: the callers pass z first, which is contiguous in ArrayNIJK
	NSUInteger w = pso.threadExecutionWidth;
	NSUInteger h = std::max<NSUInteger>(1, std::min<NSUInteger>(4, pso.maxTotalThreadsPerThreadgroup / w));
	MTLSize size = (nj>1) ? MTLSizeMake(w, h, 1) : MTLSizeMake(std::min<NSUInteger>(pso.maxTotalThreadsPerThreadgroup, 256), 1, 1);
	[enc dispatchThreads:MTLSizeMake(ni, nj, nk) threadsPerThreadgroup:size];
}

/***************************** GPU_Backend_Metal *****************************/

namespace
{
int MetalDeviceCount()
{
	@autoreleasepool
	{
		return MTLCreateSystemDefaultDevice() ? 1 : 0;
	}
}

bool MetalDeviceInfo(int n, GPU_DeviceInfo* info)
{
	@autoreleasepool
	{
		id<MTLDevice> device = MTLCreateSystemDefaultDevice();
		if ((n!=0) || !device)
			return false;
		snprintf(info->name, sizeof(info->name), "%s", [[device name] UTF8String]);
		info->memory = [device recommendedMaxWorkingSetSize];
		info->integrated = [device hasUnifiedMemory];
		return true;
	}
}

GPU_Backend* NewMetalBackend(int n)
{
	return n==0 ? GPU_Backend_Metal::New() : NULL;
}
}

const GPU_Plugin* GPU_Backend_Metal::Plugin()
{
	// built into libopenEMS, never loaded as a plugin: no build id to check
	static const GPU_Plugin plugin = {"", "metal", MetalDeviceCount, MetalDeviceInfo, NewMetalBackend};
	return &plugin;
}

GPU_Backend_Metal* GPU_Backend_Metal::New()
{
	id<MTLDevice> device = MTLCreateSystemDefaultDevice();
	if (!device)
		return NULL;
	Impl* impl = new Impl();
	impl->ctx = std::make_shared<Metal_Context>();
	impl->ctx->device = device;
	impl->ctx->queue = [device newCommandQueue];
	impl->device = device;
	impl->numCells = 0;
	return new GPU_Backend_Metal(impl);
}

GPU_Backend* GPU_Backend_Metal::NewSubGridBackend()
{
	Impl* impl = new Impl();
	impl->ctx = d->ctx;
	impl->device = d->device;
	impl->numCells = 0;
	return new GPU_Backend_Metal(impl);
}

GPU_Backend_Metal::GPU_Backend_Metal(Impl* impl)
{
	d = impl;
}

GPU_Backend_Metal::~GPU_Backend_Metal()
{
	// run out the encoded batch: its encoder must be ended, and the device must be done
	// before ARC releases the buffers; a device error was reported already, do not throw here
	try
	{
		d->Flush();
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	delete d;
}

std::string GPU_Backend_Metal::GetName() const
{
	return std::string("Metal (") + [[d->device name] UTF8String] + ")";
}

bool GPU_Backend_Metal::Init(const Operator* op, const GPU_CellModel& m)
{
	// the kernels reach the model through the gpuAddress of its buffers
	if (@available(macOS 13.0, *)) {}
	else
	{
		std::cerr << "GPU_Backend_Metal::Init: Error: the Metal backend needs macOS 13 or newer" << std::endl;
		return false;
	}
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	if ((m.numLines[0]!=numLines[0]) || (m.numLines[1]!=numLines[1]) || (m.numLines[2]!=numLines[2]))
		return false;
	d->dim.nx = numLines[0];
	d->dim.ny = numLines[1];
	d->dim.nz = numLines[2];
	d->numCells = (size_t)numLines[0]*numLines[1]*numLines[2];

	// the fields: bound to the main kernels, reached through pointers by the extensions
	const size_t bytes = 3*d->numCells*sizeof(float);
	d->volt = d->NewBuffer(bytes);
	d->curr = d->NewBuffer(bytes);
	for (id<MTLBuffer> f : {d->volt, d->curr})
	{
		d->ctx->resources.push_back(f);
		if (@available(macOS 13.0, *))
			d->ctx->addresses[f.gpuAddress] = f;
	}

	d->cm_index16 = m.Index16();
	size_t total = 0;
	for (int f=0; f<2; ++f)
	{
		const GPU_CellModel::Field& src = m.field[f];
		CM_Field& F = d->cm_field[f];
		std::memset(&F, 0, sizeof(F));
		if (d->cm_index16)
		{
			// converted straight into the shared buffer, no host copy of the whole index
			d->cm_index[f] = d->NewBuffer(src.index.size()*sizeof(uint16_t));
			std::copy(src.index.begin(), src.index.end(), static_cast<uint16_t*>([d->cm_index[f] contents]));
		}
		else
			d->cm_index[f] = d->NewBuffer(src.index.size()*sizeof(uint32_t), src.index.data());
		auto res = [&](const auto& v) {total += v.size()*sizeof(v[0]); return d->NewResource(v.size()*sizeof(v[0]), v.empty() ? NULL : v.data());};
		F.combos = (const CM_Combo*)res(src.combos);
		F.slot = (const uint32_t*)res(src.slot);
		F.numSlots = src.numSlots;
		F.state = (float*)d->NewResource(src.numSlots*src.slotFloats*sizeof(float));   // zero
		total += src.numSlots*src.slotFloats*sizeof(float);
		F.poles = (const CM_Pole*)res(src.poles);
		for (int q=0; q<5; ++q)
			F.profile[q] = (const float*)res(src.profile[q]);
		for (int n=0; n<3; ++n)
		{
			F.gxy[n] = (const float*)res(src.gxy[n]);
			F.gz[n] = (const float*)res(src.gz[n]);
		}
		if (g_settings.GetVerboseLevel()>0)
			std::cout << "GPU_Backend_Metal: cell model, " << (f ? "currents" : "voltages") << ": " << src.combos.size() << " combinations, "
			          << src.numSlots << " state slots of " << src.slotFloats << " floats, " << src.poles.size() << " poles" << std::endl;
	}
	if (g_settings.GetVerboseLevel()>0)
		std::cout << "GPU_Backend_Metal: cell model on the device, " << (total + 2*d->numCells*(d->cm_index16 ? 2 : 4))/1048576
		          << " MiB, " << (d->cm_index16 ? 16 : 32) << " bit node index" << std::endl;

	// compile the kernels now, not in the first timestep
	d->Pipeline(d->cm_index16 ? "cm_voltages16" : "cm_voltages32");
	d->Pipeline(d->cm_index16 ? "cm_currents16" : "cm_currents32");
	return true;
}

void GPU_Backend_Metal::Update(bool voltages, unsigned int x0, unsigned int x1)
{
	x1 = std::min(x1, d->dim.nx);
	if (x0>=x1)
		return;
	const char* name = voltages ? (d->cm_index16 ? "cm_voltages16" : "cm_voltages32") : (d->cm_index16 ? "cm_currents16" : "cm_currents32");
	id<MTLComputePipelineState> pso = d->Pipeline(name);
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:(voltages ? d->volt : d->curr) offset:0 atIndex:0];
	[enc setBuffer:(voltages ? d->curr : d->volt) offset:0 atIndex:1];
	[enc setBuffer:d->cm_index[voltages ? 0 : 1] offset:0 atIndex:2];
	[enc setBytes:&d->cm_field[voltages ? 0 : 1] length:sizeof(CM_Field) atIndex:3];
	const uint32_t grid[4] = {d->dim.nx, d->dim.ny, d->dim.nz, x0};   // CM_Grid
	[enc setBytes:grid length:sizeof(grid) atIndex:4];
	d->Dispatch(pso, d->dim.nz, d->dim.ny, x1-x0);
}

void GPU_Backend_Metal::UpdateVoltages(const Engine_Box& box) {Update(true, box.start[0], box.stop[0]);}
void GPU_Backend_Metal::UpdateCurrents(const Engine_Box& box) {Update(false, box.start[0], box.stop[0]);}

void* GPU_Backend_Metal::DeviceAlloc(size_t bytes, const void* host)
{
	return (void*)d->NewResource(bytes, host);
}

//! the buffer of DeviceAlloc() that holds the device address \a dev, and the offset in it
static id<MTLBuffer> BufferOf(const std::map<uint64_t, id<MTLBuffer>>& addresses, const void* dev, size_t& offset)
{
	const uint64_t a = (uint64_t)dev;
	std::map<uint64_t, id<MTLBuffer>>::const_iterator it = addresses.upper_bound(a);
	if (it==addresses.begin())
		throw std::logic_error("GPU_Backend_Metal: not a device address of the backend");
	--it;
	offset = a - it->first;
	if (offset>=[it->second length])
		throw std::logic_error("GPU_Backend_Metal: not a device address of the backend");
	return it->second;
}

void GPU_Backend_Metal::DeviceDownload(void* host, const void* dev, size_t bytes)
{
	size_t offset;
	id<MTLBuffer> buf = BufferOf(d->ctx->addresses, dev, offset);
	d->Flush();
	std::memcpy(host, (const char*)[buf contents] + offset, bytes);
}

float* GPU_Backend_Metal::DeviceField(bool currents)
{
	uint64_t address = 0;
	if (@available(macOS 13.0, *))
		address = (currents ? d->curr : d->volt).gpuAddress;
	return (float*)address;
}

void GPU_Backend_Metal::RunPoints(GPU_PointOp kind, const void* op, uint64_t begin, uint64_t end, bool serial)
{
	if (end<=begin)
		return;
	const char* name = NULL;
	size_t size = 0;
	switch (kind)
	{
#define OP(K, T, n) case K: name = serial ? #n "_serial" : #n; size = sizeof(T); break;
	CM_POINT_OPS(OP)
#undef OP
	}
	if (name==NULL)
		throw std::logic_error("GPU_Backend_Metal: unknown point operation");
	id<MTLComputePipelineState> pso = d->Pipeline(name);
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBytes:op length:size atIndex:0];
	const uint64_t range[2] = {begin, end};   // ulong2
	[enc setBytes:range length:sizeof(range) atIndex:1];
	d->Dispatch(pso, serial ? 1 : end-begin);
}

static void CopyField(const void* src, void* dst, size_t src_bytes, size_t dst_bytes)
{
	if (src_bytes!=dst_bytes)
		throw std::runtime_error("GPU_Backend_Metal: host and device field size mismatch");
	if (src!=dst)   // nothing to do if the host mirror is the shared buffer
		std::memcpy(dst, src, src_bytes);
}

// unified memory: Engine_GPU maps its host mirror onto these buffers, so the fields never
// travel and Download() only has to wait for the device
FDTD_FLOAT* GPU_Backend_Metal::GetShared(bool currents) const
{
	return static_cast<FDTD_FLOAT*>([(currents ? d->curr : d->volt) contents]);
}

// the host reads the shared buffers directly, so it has to wait for the encoded batch
void GPU_Backend_Metal::Synchronize()
{
	d->Flush();
}

// The dump thread reads the fields (and sums the frequency domain dumps) while the engine keeps
// stepping, so the device copies them into buffers of their own instead of sharing the field buffers.
bool GPU_Backend_Metal::SnapshotFields(const FDTD_FLOAT* &volt, const FDTD_FLOAT* &curr)
{
	if (d->snap_failed)
		return false;
	const size_t bytes = 3*d->numCells*sizeof(float);
	if (!d->snap_volt)
	{
		// optional, allocated during the run: without them the dumps read the shared field buffers
		d->snap_volt = [d->device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
		d->snap_curr = [d->device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
		if (!d->snap_volt || !d->snap_curr)
		{
			d->snap_volt = nil;
			d->snap_curr = nil;
			d->snap_failed = true;
			std::cout << "GPU_Backend_Metal: not enough memory for field snapshots, the dumps read the fields" << std::endl;
			return false;
		}
	}
	id<MTLComputePipelineState> pso = d->Pipeline("copy_fields");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:d->curr offset:0 atIndex:1];
	[enc setBuffer:d->snap_volt offset:0 atIndex:2];
	[enc setBuffer:d->snap_curr offset:0 atIndex:3];
	const uint32_t grid[4] = {d->dim.nx, d->dim.ny, d->dim.nz, 0};   // CM_Grid
	[enc setBytes:grid length:sizeof(grid) atIndex:4];
	d->Dispatch(pso, d->dim.nz, d->dim.ny, 3*(size_t)d->dim.nx);
	// commit without waiting: the next timesteps are encoded into a new command buffer,
	// which the device executes after this one
	d->snap_cmd = d->ctx->Commit();
	volt = static_cast<const FDTD_FLOAT*>([d->snap_volt contents]);
	curr = static_cast<const FDTD_FLOAT*>([d->snap_curr contents]);
	return true;
}

// called from the dump thread, while the main thread encodes the next timesteps
void GPU_Backend_Metal::WaitSnapshot()
{
	id<MTLCommandBuffer> cmd = d->snap_cmd;
	if (cmd)
		[cmd waitUntilCompleted];
}

bool GPU_Backend_Metal::CalcFastEnergy(const unsigned int numNodes[3], double& E_energy, double& H_energy)
{
	E_energy = H_energy = 0;
	Metal_GridDim L = {numNodes[0], numNodes[1], numNodes[2]};
	if ((L.nx>d->dim.nx) || (L.ny>d->dim.ny) || (L.nz>d->dim.nz))
		return false;
	const size_t count = (size_t)L.ny*L.nz;
	if (L.nx==0 || count==0)
		return true;
	if (!d->energy || ([d->energy length] < count*2*sizeof(float)))
		d->energy = d->NewBuffer(count*2*sizeof(float));

	id<MTLComputePipelineState> pso = d->Pipeline("energy_lines");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:d->curr offset:0 atIndex:1];
	[enc setBuffer:d->energy offset:0 atIndex:2];
	const uint32_t grid[4] = {d->dim.nx, d->dim.ny, d->dim.nz, 0}, lines[4] = {L.nx, L.ny, L.nz, 0};   // CM_Grid
	[enc setBytes:grid length:sizeof(grid) atIndex:3];
	[enc setBytes:lines length:sizeof(lines) atIndex:4];
	d->Dispatch(pso, L.nz, L.ny);
	d->Flush();

	const float* partial = static_cast<const float*>([d->energy contents]);
	for (size_t k=0; k<count; ++k)
	{
		E_energy += partial[2*k];
		H_energy += partial[2*k+1];
	}
	return true;
}

void GPU_Backend_Metal::Download(bool currents, ArrayLib::ArrayNIJK<FDTD_FLOAT>& field)
{
	d->Flush();
	CopyField([(currents ? d->curr : d->volt) contents], field.data(), 3*d->numCells*sizeof(float), field.size()*sizeof(FDTD_FLOAT));
}
