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
#include <cstring>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "metal_internal.h"
#include "FDTD/operator.h"

// Main FDTD updates, operation by operation the same as Engine::UpdateVoltages/UpdateCurrents.
// One thread per mesh node, gid = (z, y, x) relative to the start node O of the update range.
static const char* BASE_SOURCE = R"MSL(
kernel void update_voltages(device float* volt       [[buffer(0)]],
                            const device float* curr [[buffer(1)]],
                            const device float* vv   [[buffer(2)]],
                            const device float* vi   [[buffer(3)]],
                            constant GridDim& N      [[buffer(4)]],
                            constant GridDim& O      [[buffer(5)]],
                            uint3 gid [[thread_position_in_grid]])
{
	const uint z = O.nz + gid.x, y = O.ny + gid.y, x = O.nx + gid.z;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const uint xm = (x>0) ? N.ny*N.nz : 0;   // shift to the previous line, none on the first
	const uint ym = (y>0) ? N.nz : 0;
	const uint zm = (z>0) ? 1 : 0;
	float v;

	//for x
	v  = volt[i] * vv[i];
	v += vi[i] * (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	volt[i] = v;

	//for y
	v  = volt[sn+i] * vv[sn+i];
	v += vi[sn+i] * (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	volt[sn+i] = v;

	//for z
	v  = volt[2*sn+i] * vv[2*sn+i];
	v += vi[2*sn+i] * (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	volt[2*sn+i] = v;
}

kernel void update_currents(device float* curr       [[buffer(0)]],
                            const device float* volt [[buffer(1)]],
                            const device float* ii   [[buffer(2)]],
                            const device float* iv   [[buffer(3)]],
                            constant GridDim& N      [[buffer(4)]],
                            constant GridDim& O      [[buffer(5)]],
                            uint3 gid [[thread_position_in_grid]])
{
	const uint z = O.nz + gid.x, y = O.ny + gid.y, x = O.nx + gid.z;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const uint xp = N.ny*N.nz;
	const uint yp = N.nz;
	float c;

	//for x
	c  = curr[i] * ii[i];
	c += iv[i] * (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
	curr[i] = c;

	//for y
	c  = curr[sn+i] * ii[sn+i];
	c += iv[sn+i] * (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
	curr[sn+i] = c;

	//for z
	c  = curr[2*sn+i] * ii[2*sn+i];
	c += iv[2*sn+i] * (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	curr[2*sn+i] = c;
}

// Same with compressed coefficients: the 12 coefficients of a node are the set
// coeff[12*index[node]], in the order vv[3], vi[3], ii[3], iv[3].
kernel void update_voltages_c(device float* volt        [[buffer(0)]],
                              const device float* curr  [[buffer(1)]],
                              const device ushort* index[[buffer(2)]],
                              const device float* coeff [[buffer(3)]],
                              constant GridDim& N       [[buffer(4)]],
                              constant GridDim& O       [[buffer(5)]],
                              uint3 gid [[thread_position_in_grid]])
{
	const uint z = O.nz + gid.x, y = O.ny + gid.y, x = O.nx + gid.z;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const uint xm = (x>0) ? N.ny*N.nz : 0;
	const uint ym = (y>0) ? N.nz : 0;
	const uint zm = (z>0) ? 1 : 0;
	const device float* vv = coeff + 12*uint(index[i]);
	const device float* vi = vv + 3;
	float v;

	//for x
	v  = volt[i] * vv[0];
	v += vi[0] * (curr[2*sn+i] - curr[2*sn+i-ym] - curr[sn+i] + curr[sn+i-zm]);
	volt[i] = v;

	//for y
	v  = volt[sn+i] * vv[1];
	v += vi[1] * (curr[i] - curr[i-zm] - curr[2*sn+i] + curr[2*sn+i-xm]);
	volt[sn+i] = v;

	//for z
	v  = volt[2*sn+i] * vv[2];
	v += vi[2] * (curr[sn+i] - curr[sn+i-xm] - curr[i] + curr[i-ym]);
	volt[2*sn+i] = v;
}

kernel void update_currents_c(device float* curr        [[buffer(0)]],
                              const device float* volt  [[buffer(1)]],
                              const device ushort* index[[buffer(2)]],
                              const device float* coeff [[buffer(3)]],
                              constant GridDim& N       [[buffer(4)]],
                              constant GridDim& O       [[buffer(5)]],
                              uint3 gid [[thread_position_in_grid]])
{
	const uint z = O.nz + gid.x, y = O.ny + gid.y, x = O.nx + gid.z;
	const uint sn = N.nx*N.ny*N.nz;
	const uint i  = nijk(N, 0, x, y, z);
	const uint xp = N.ny*N.nz;
	const uint yp = N.nz;
	const device float* ii = coeff + 12*uint(index[i]) + 6;
	const device float* iv = ii + 3;
	float c;

	//for x
	c  = curr[i] * ii[0];
	c += iv[0] * (volt[2*sn+i] - volt[2*sn+i+yp] - volt[sn+i] + volt[sn+i+1]);
	curr[i] = c;

	//for y
	c  = curr[sn+i] * ii[1];
	c += iv[1] * (volt[i] - volt[i+1] - volt[2*sn+i] + volt[2*sn+i+xp]);
	curr[sn+i] = c;

	//for z
	c  = curr[2*sn+i] * ii[2];
	c += iv[2] * (volt[sn+i] - volt[sn+i+xp] - volt[i] + volt[i+yp]);
	curr[2*sn+i] = c;
}

// Field energy: the squared voltages and currents of the nodes below L, summed along x
// in float, one thread per (y,z) line; gid = (z, y). The host sums the lines in double.
kernel void energy_lines(const device float* volt [[buffer(0)]],
                         const device float* curr [[buffer(1)]],
                         device float2* partial   [[buffer(2)]],
                         constant GridDim& N      [[buffer(3)]],
                         constant GridDim& L      [[buffer(4)]],
                         uint2 gid [[thread_position_in_grid]])
{
	const uint z = gid.x, y = gid.y;
	const uint sn = N.nx*N.ny*N.nz;
	float E = 0, H = 0;
	for (uint x=0; x<L.nx; ++x)
	{
		const uint i = nijk(N, 0, x, y, z);
		E += volt[i]*volt[i];
		E += volt[sn+i]*volt[sn+i];
		E += volt[2*sn+i]*volt[2*sn+i];
		H += curr[i]*curr[i];
		H += curr[sn+i]*curr[sn+i];
		H += curr[2*sn+i]*curr[2*sn+i];
	}
	partial[y*L.nz + z] = float2(E, H);
}
)MSL";

/***************************** Impl *****************************/

id<MTLComputePipelineState> Metal_Context::Pipeline(const char* source, const char* function)
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

		NSString* src = [NSString stringWithFormat:@"%s\n%s", METAL_COMMON_SOURCE, source];
		NSError* error = nil;
		id<MTLLibrary> lib = [device newLibraryWithSource:src options:options error:&error];
		if (!lib)
			throw std::runtime_error(std::string("GPU_Backend_Metal: kernel compilation failed: ") + [[error localizedDescription] UTF8String]);
		id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:function]];
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
	@autoreleasepool
	{
		cmd = [queue commandBuffer];
		// concurrent dispatch with explicit barriers, see GPU_Backend_Metal::Impl::Dispatch()
		enc = [cmd computeCommandEncoderWithDispatchType:MTLDispatchTypeConcurrent];
		group = NULL;
	}
	return enc;
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
		if (cmd.status==MTLCommandBufferStatusError)
			throw std::runtime_error(std::string("GPU_Backend_Metal: command buffer failed: ") + [[cmd.error localizedDescription] UTF8String]);
		enc = nil;
		cmd = nil;
	}
}

id<MTLBuffer> GPU_Backend_Metal::Impl::NewBuffer(size_t bytes, const void* data)
{
	id<MTLBuffer> buf = [device newBufferWithLength:std::max(bytes, (size_t)4) options:MTLResourceStorageModeShared];
	if (!buf)
		throw std::runtime_error("GPU_Backend_Metal: buffer allocation failed");
	std::memset([buf contents], 0, [buf length]);
	if (data && bytes)
		std::memcpy([buf contents], data, bytes);
	return buf;
}

void GPU_Backend_Metal::Impl::Dispatch(id<MTLComputePipelineState> pso, size_t ni, size_t nj, size_t nk, const void* group)
{
	if (ni==0 || nj==0 || nk==0)
		return;
	id<MTLComputeCommandEncoder> enc = Encoder();
	if ((group==NULL) || (group!=ctx->group))
		[enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
	ctx->group = group;
	NSUInteger w = pso.threadExecutionWidth;
	NSUInteger h = std::max<NSUInteger>(1, std::min<NSUInteger>(4, pso.maxTotalThreadsPerThreadgroup / w));
	MTLSize size = (nj>1) ? MTLSizeMake(w, h, 1) : MTLSizeMake(std::min<NSUInteger>(pso.maxTotalThreadsPerThreadgroup, 256), 1, 1);
	[enc dispatchThreads:MTLSizeMake(ni, nj, nk) threadsPerThreadgroup:size];
}

void GPU_Backend_Metal::Impl::SetGridDim(unsigned int index)
{
	[Encoder() setBytes:&dim length:sizeof(dim) atIndex:index];
}

/***************************** GPU_Backend_Metal *****************************/

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
	d->Flush();
	delete d;
}

std::string GPU_Backend_Metal::GetName() const
{
	return std::string("Metal (") + [[d->device name] UTF8String] + ")";
}

bool GPU_Backend_Metal::Init(const Operator* op)
{
	unsigned int numLines[3];
	for (int n=0; n<3; ++n)
		numLines[n] = op->GetNumberOfLines(n, true);
	d->dim.nx = numLines[0];
	d->dim.ny = numLines[1];
	d->dim.nz = numLines[2];
	d->numCells = (size_t)numLines[0]*numLines[1]*numLines[2];
	d->main_start.nx = d->main_start.ny = d->main_start.nz = 0;
	d->main_stop = d->dim;
	d->upml_fused = -1;

	// the kernels index with 32 bit
	if (3*d->numCells > std::numeric_limits<uint32_t>::max())
	{
		std::cerr << "GPU_Backend_Metal::Init: Error: the mesh is too large for 32 bit indexing" << std::endl;
		return false;
	}

	const size_t bytes = 3*d->numCells*sizeof(float);
	d->volt = d->NewBuffer(bytes);
	d->curr = d->NewBuffer(bytes);

	// upload the final operator coefficients, including all changes by operator extensions
	std::vector<float> vv(3*d->numCells), vi(3*d->numCells), ii(3*d->numCells), iv(3*d->numCells);
	unsigned int pos[3];
	size_t idx = 0;
	for (int n=0; n<3; ++n)
		for (pos[0]=0; pos[0]<numLines[0]; ++pos[0])
			for (pos[1]=0; pos[1]<numLines[1]; ++pos[1])
				for (pos[2]=0; pos[2]<numLines[2]; ++pos[2], ++idx)
				{
					vv[idx] = op->GetVV(n, pos[0], pos[1], pos[2]);
					vi[idx] = op->GetVI(n, pos[0], pos[1], pos[2]);
					ii[idx] = op->GetII(n, pos[0], pos[1], pos[2]);
					iv[idx] = op->GetIV(n, pos[0], pos[1], pos[2]);
				}

	if (CompressCoefficients(vv, vi, ii, iv))
	{
		d->Pipeline(BASE_SOURCE, "update_voltages_c");
		d->Pipeline(BASE_SOURCE, "update_currents_c");
		return true;
	}

	d->vv = d->NewBuffer(bytes, vv.data());
	d->vi = d->NewBuffer(bytes, vi.data());
	d->ii = d->NewBuffer(bytes, ii.data());
	d->iv = d->NewBuffer(bytes, iv.data());

	// compile the main kernels now, not in the first timestep
	d->Pipeline(BASE_SOURCE, "update_voltages");
	d->Pipeline(BASE_SOURCE, "update_currents");
	return true;
}

namespace
{
//! bit pattern of the 12 coefficients of a node, see update_voltages_c
struct CoeffSet
{
	uint32_t v[12];
	bool operator==(const CoeffSet& o) const {return std::memcmp(v, o.v, sizeof(v))==0;}
};

struct CoeffSetHash
{
	size_t operator()(const CoeffSet& s) const
	{
		uint64_t h = 1469598103934665603ULL;   // FNV-1a
		for (int k=0; k<12; ++k)
			h = (h ^ s.v[k]) * 1099511628211ULL;
		return h;
	}
};
}

// Most nodes share one of a few coefficient sets (same material and mesh spacing).
// Store every distinct set once and a 16 bit set index per node, which removes most
// of the coefficient memory traffic of the main updates. Bit-exact, since the sets
// hold the original values. Not used if there are more sets than a 16 bit index holds.
bool GPU_Backend_Metal::CompressCoefficients(const std::vector<float>& vv, const std::vector<float>& vi,
                                             const std::vector<float>& ii, const std::vector<float>& iv)
{
	const size_t sn = d->numCells;
	const size_t max_sets = (size_t)std::numeric_limits<uint16_t>::max()+1;
	std::vector<uint16_t> index(sn);
	std::vector<float> table;
	std::unordered_map<CoeffSet, uint16_t, CoeffSetHash> sets;

	const std::vector<float>* src[4] = {&vv, &vi, &ii, &iv};
	CoeffSet prev;
	bool have_prev = false;
	for (size_t i=0; i<sn; ++i)
	{
		CoeffSet s;
		for (int c=0; c<4; ++c)
			for (int n=0; n<3; ++n)
				std::memcpy(&s.v[3*c+n], &(*src[c])[n*sn+i], sizeof(float));
		// neighbours along z mostly share the set
		if (have_prev && (s==prev))
		{
			index[i] = index[i-1];
			continue;
		}
		std::unordered_map<CoeffSet, uint16_t, CoeffSetHash>::const_iterator it = sets.find(s);
		if (it==sets.end())
		{
			if (sets.size()>=max_sets)
				return false;
			it = sets.insert(std::make_pair(s, (uint16_t)sets.size())).first;
			table.insert(table.end(), (const float*)s.v, (const float*)s.v+12);
		}
		index[i] = it->second;
		prev = s;
		have_prev = true;
	}

	d->index = d->NewBuffer(sn*sizeof(uint16_t), index.data());
	d->coeff = d->NewBuffer(table.size()*sizeof(float), table.data());
	std::cout << "GPU_Backend_Metal: " << sets.size() << " distinct coefficient sets, compressed update coefficients" << std::endl;
	return true;
}

void GPU_Backend_Metal::UpdateVoltages()
{
	id<MTLComputePipelineState> pso = d->Pipeline(BASE_SOURCE, d->coeff ? "update_voltages_c" : "update_voltages");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:d->curr offset:0 atIndex:1];
	[enc setBuffer:(d->coeff ? d->index : d->vv) offset:0 atIndex:2];
	[enc setBuffer:(d->coeff ? d->coeff : d->vi) offset:0 atIndex:3];
	d->SetGridDim(4);
	[enc setBytes:&d->main_start length:sizeof(d->main_start) atIndex:5];
	const Metal_GridDim& b = d->main_start;
	const Metal_GridDim& e = d->main_stop;
	d->Dispatch(pso, e.nz-b.nz, e.ny-b.ny, e.nx-b.nx);
}

void GPU_Backend_Metal::UpdateCurrents()
{
	id<MTLComputePipelineState> pso = d->Pipeline(BASE_SOURCE, d->coeff ? "update_currents_c" : "update_currents");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->curr offset:0 atIndex:0];
	[enc setBuffer:d->volt offset:0 atIndex:1];
	[enc setBuffer:(d->coeff ? d->index : d->ii) offset:0 atIndex:2];
	[enc setBuffer:(d->coeff ? d->coeff : d->iv) offset:0 atIndex:3];
	d->SetGridDim(4);
	[enc setBytes:&d->main_start length:sizeof(d->main_start) atIndex:5];
	// the currents on the last mesh line are not updated
	const Metal_GridDim& b = d->main_start;
	const unsigned int ex = std::min(d->main_stop.nx, d->dim.nx-1);
	const unsigned int ey = std::min(d->main_stop.ny, d->dim.ny-1);
	const unsigned int ez = std::min(d->main_stop.nz, d->dim.nz-1);
	if ((ex>b.nx) && (ey>b.ny) && (ez>b.nz))
		d->Dispatch(pso, ez-b.nz, ey-b.ny, ex-b.nx);
}

static void CopyField(const void* src, void* dst, size_t src_bytes, size_t dst_bytes)
{
	if (src_bytes!=dst_bytes)
		throw std::runtime_error("GPU_Backend_Metal: host and device field size mismatch");
	if (src!=dst)   // nothing to do if the host mirror is the shared buffer
		std::memcpy(dst, src, src_bytes);
}

FDTD_FLOAT* GPU_Backend_Metal::GetSharedVoltages() const
{
	return static_cast<FDTD_FLOAT*>([d->volt contents]);
}

FDTD_FLOAT* GPU_Backend_Metal::GetSharedCurrents() const
{
	return static_cast<FDTD_FLOAT*>([d->curr contents]);
}

void GPU_Backend_Metal::Synchronize()
{
	d->Flush();
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

	id<MTLComputePipelineState> pso = d->Pipeline(BASE_SOURCE, "energy_lines");
	id<MTLComputeCommandEncoder> enc = d->Encoder();
	[enc setComputePipelineState:pso];
	[enc setBuffer:d->volt offset:0 atIndex:0];
	[enc setBuffer:d->curr offset:0 atIndex:1];
	[enc setBuffer:d->energy offset:0 atIndex:2];
	d->SetGridDim(3);
	[enc setBytes:&L length:sizeof(L) atIndex:4];
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

void GPU_Backend_Metal::DownloadVoltages(ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	d->Flush();
	CopyField([d->volt contents], volt.data(), 3*d->numCells*sizeof(float), volt.size()*sizeof(FDTD_FLOAT));
}

void GPU_Backend_Metal::DownloadCurrents(ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	d->Flush();
	CopyField([d->curr contents], curr.data(), 3*d->numCells*sizeof(float), curr.size()*sizeof(FDTD_FLOAT));
}

void GPU_Backend_Metal::UploadVoltages(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& volt)
{
	d->Flush();
	CopyField(volt.data(), [d->volt contents], volt.size()*sizeof(FDTD_FLOAT), 3*d->numCells*sizeof(float));
}

void GPU_Backend_Metal::UploadCurrents(const ArrayLib::ArrayNIJK<FDTD_FLOAT>& curr)
{
	d->Flush();
	CopyField(curr.data(), [d->curr contents], curr.size()*sizeof(FDTD_FLOAT), 3*d->numCells*sizeof(float));
}

// all Metal extensions, see metal_internal.h
static const Metal_ExtensionFactory METAL_EXTENSIONS[] = {
	Metal_CreateExt_Excitation,
	Metal_CreateExt_UPML,
	Metal_CreateExt_Mur_ABC,
	Metal_CreateExt_LorentzMaterial,
	Metal_CreateExt_LumpedRLC,
	Metal_CreateExt_TFSF,
	Metal_CreateExt_Absorbing_BC,
	Metal_CreateExt_SteadyState,
	Metal_CreateExt_Cylinder,
};

GPU_Extension* GPU_Backend_Metal::CreateExtension(Engine_Extension* eng_ext, Engine* eng)
{
	for (size_t n=0; n<sizeof(METAL_EXTENSIONS)/sizeof(METAL_EXTENSIONS[0]); ++n)
	{
		GPU_Extension* gpu_ext = METAL_EXTENSIONS[n](d, eng_ext, eng);
		if (gpu_ext)
			return gpu_ext;
	}
	return NULL;
}
