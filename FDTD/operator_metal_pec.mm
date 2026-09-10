/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "operator_metal.h"
#include "ContinuousStructure.h"
#include "CSPrimPolygon.h"
#include "CSPrimLinPoly.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
// Exact FP64 Yee-index bounds plus conservative FP32 polygon predicates.
// Near edges or unsupported primitives, resolve the query with CSXCAD FP64.
struct PecPrimitive
{
	uint32_t bounds[12]; // [axis][primal/dual][begin/end), computed in FP64
	uint32_t kind, normal, first, count;
};
struct PecPoint { float x, y; };
const uint32_t cpuQuery = UINT32_MAX - 1;
const char* pecSource = R"METAL(
#include <metal_stdlib>
using namespace metal;
struct Primitive { uint bounds[12]; uint kind, normal, first, count; };
kernel void pec_mask(const device Primitive* prims [[buffer(0)]],
                     const device float2* vertices [[buffer(1)]],
                     const device uint* offsets [[buffer(2)]],
                     const device uint* candidates [[buffer(3)]],
                     const device float* grid [[buffer(4)]],
                     device uint* winners [[buffer(5)]],
                     constant uint4& p [[buffer(6)]],
                     uint3 gid [[thread_position_in_grid]])
{
	if (gid.x >= p.z * 3 || gid.y >= p.y) return;
	uint n = gid.x % 3, z = gid.x / 3, y = gid.y, x = p.w;
	uint line = y;
	uint output = line * p.z * 3 + gid.x;
	float c[3] = {grid[2*x + (n == 0)],
	              grid[2*p.x + 2*y + (n == 1)],
	              grid[2*(p.x+p.y) + 2*z + (n == 2)]};
	uint result = 0xffffffffu;
	for (uint k = offsets[line]; k < offsets[line+1]; ++k)
	{
		uint id = candidates[k];
		Primitive q = prims[id];
		if (q.kind == 0) { result = 0xfffffffeu; break; }
		bool outside = false, uncertain = false;
		uint pos[3] = {x,y,z};
		for (uint a = 0; a < 3; ++a)
		{
			uint b = a*4 + 2*(a == n);
			outside |= pos[a] < q.bounds[b] || pos[a] >= q.bounds[b+1];
		}
		if (outside) continue;
		if (q.kind == 1) { result = id; break; }
		float px = c[(q.normal+1)%3], py = c[(q.normal+2)%3];
		float2 prev = vertices[q.first+q.count-1];
		int winding = 0;
		bool over = prev.y >= py;
		for (uint j = 0; j < q.count; ++j)
		{
			float2 v = vertices[q.first+j];
			float scale = max(1.0f, max(max(abs(px),abs(py)), max(max(abs(v.x),abs(v.y)), max(abs(prev.x),abs(prev.y)))));
			float e = 32.0f * FLT_EPSILON * scale;
			if (abs(v.y-py) <= e || abs(prev.y-py) <= e) { uncertain = true; break; }
			bool nextOver = v.y >= py;
			if (over != nextOver)
			{
				float left = (v.y-py)*(v.x-prev.x);
				float right = (v.y-prev.y)*(v.x-px);
				if (abs(left-right) <= 128.0f*FLT_EPSILON*scale*scale) { uncertain = true; break; }
				if (left <= right) { if (nextOver) ++winding; }
				else if (!nextOver) --winding;
			}
			over = nextOver; prev = v;
		}
		if (uncertain) { result = 0xfffffffeu; break; }
		if (winding != 0) { result = id; break; }
	}
	winners[output] = result;
}
)METAL";
}

bool Operator_Metal::CalcPEC()
{
	const auto started = std::chrono::steady_clock::now();
	const char* setting = std::getenv("OPENEMS_METAL_PEC");
	// --engine=metal enables PEC mapping too; retain a diagnostic CPU override.
	if ((setting && setting[0] == '0') || m_MeshType != CARTESIAN)
	{
		const bool result = Operator::CalcPEC();
		if (setting)
			std::cout << "CPU PEC: " << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << " s" << std::endl;
		return result;
	}
	const bool verify = setting && std::string(setting) == "verify";
	@autoreleasepool
	{
		id<MTLDevice> device = MTLCreateSystemDefaultDevice();
		id<MTLCommandQueue> queue = [device newCommandQueue];
		NSError* error = nil;
		MTLCompileOptions* options = [MTLCompileOptions new];
		options.fastMathEnabled = NO; // Geometry decisions must not use field fast math.
		id<MTLLibrary> library = [device newLibraryWithSource:@(pecSource) options:options error:&error];
		id<MTLFunction> function = [library newFunctionWithName:@"pec_mask"];
		id<MTLComputePipelineState> pipeline = function ? [device newComputePipelineStateWithFunction:function error:&error] : nil;
		if (!queue || !pipeline)
		{
			std::cerr << "Metal PEC: unavailable; using CPU: " << (error ? error.localizedDescription.UTF8String : "no device") << std::endl;
			return Operator::CalcPEC();
		}
		auto buffer = [device](const void* data, size_t bytes) {
			id<MTLBuffer> b = data && bytes ? [device newBufferWithBytes:data length:bytes options:MTLResourceStorageModeShared]
			    : [device newBufferWithLength:std::max<size_t>(bytes, sizeof(PecPrimitive)) options:MTLResourceStorageModeShared];
			if (!b) throw std::runtime_error("Metal PEC: buffer allocation failed");
			return b;
		};
		const auto types = static_cast<CSProperties::PropertyType>(CSProperties::MATERIAL | CSProperties::METAL);
		const auto primitives = CSX->GetAllPrimitives(true, types);
		if (primitives.size() >= cpuQuery) return Operator::CalcPEC();
		// Convert bounds to exact Yee-index ranges on the CPU. In particular,
		// zero-thickness sheets must not become a tolerance-thickened volume.
		std::vector<double> axes[3][2];
		std::vector<float> grid;
		for (int a = 0; a < 3; ++a)
			for (unsigned i = 0; i < numLines[a]; ++i)
				for (int dual = 0; dual < 2; ++dual)
				{
					double coord = GetDiscLine(a, i, dual);
					if (!std::isfinite(coord) || std::abs(coord) > 1e10) return Operator::CalcPEC();
					axes[a][dual].push_back(coord);
					grid.push_back(static_cast<float>(coord));
				}
		for (int a = 0; a < 3; ++a)
			for (int dual = 0; dual < 2; ++dual)
				if (!std::is_sorted(axes[a][dual].begin(), axes[a][dual].end())) return Operator::CalcPEC();
		std::vector<PecPrimitive> flattened;
		std::vector<PecPoint> vertices;
		size_t unsupported = 0;
		for (auto* prim : primitives)
		{
			PecPrimitive q{};
			const int type = prim->GetType();
			if (!prim->HasTransform() && prim->GetCoordInputType() == CARTESIAN &&
			    (prim->GetCoordinateSystem() == CARTESIAN || prim->GetCoordinateSystem() == UNDEFINED_CS) &&
			    (type == CSPrimitives::BOX || type == CSPrimitives::POLYGON || type == CSPrimitives::LINPOLY))
			{
				double bounds[6];
				prim->GetBoundBox(bounds);
				q.kind = type == CSPrimitives::BOX ? 1 : 2;
				for (int a = 0; a < 3; ++a)
				{
					if (!std::isfinite(bounds[2*a]) || !std::isfinite(bounds[2*a+1])) q.kind = 0;
					for (int dual = 0; dual < 2; ++dual)
					{
						const auto& axis = axes[a][dual];
						q.bounds[a*4+dual*2] = std::lower_bound(axis.begin(), axis.end(), bounds[2*a])-axis.begin();
						q.bounds[a*4+dual*2+1] = std::upper_bound(axis.begin(), axis.end(), bounds[2*a+1])-axis.begin();
					}
				}
				if (q.kind == 2)
				{
					auto* polygon = static_cast<CSPrimPolygon*>(prim);
					if (polygon->GetQtyCoords() > UINT32_MAX - vertices.size())
						return Operator::CalcPEC();
					q.normal = polygon->GetNormDir();
					q.first = static_cast<uint32_t>(vertices.size());
					q.count = static_cast<uint32_t>(polygon->GetQtyCoords());
					if (!q.count || q.normal > 2) q.kind = 0;
					for (uint32_t j = 0; j < q.count; ++j)
					{
						double vx = polygon->GetCoord(2*j), vy = polygon->GetCoord(2*j+1);
						if (!std::isfinite(vx) || !std::isfinite(vy) || std::abs(vx) > 1e10 || std::abs(vy) > 1e10) q.kind = 0;
						vertices.push_back({static_cast<float>(vx), static_cast<float>(vy)});
					}
				}
			}
			if (!q.kind) ++unsupported;
			flattened.push_back(q);
		}
		if (unsupported)
			std::cout << "Metal PEC: " << unsupported << " unsupported primitives; affected queries use CPU" << std::endl;
		id<MTLBuffer> geometry = buffer(flattened.data(), flattened.size()*sizeof(PecPrimitive));
		id<MTLBuffer> points = buffer(vertices.data(), vertices.size()*sizeof(PecPoint));
		id<MTLBuffer> coordinates = buffer(grid.data(), grid.size()*sizeof(float));
		size_t resolved = 0, queries = 0;
		std::fill(m_Nr_PEC, m_Nr_PEC+3, 0);
		// One X slab bounds mask/index memory. Keep CSXCAD's exact candidate order,
		// including equal priorities and its boundary-domain culling behavior.
		for (unsigned x = 0; x < numLines[0]; ++x)
		@autoreleasepool
		{
			std::vector<uint32_t> slab;
			for (uint32_t id = 0; id < flattened.size(); ++id)
			{
				const auto& q = flattened[id];
				if (!q.kind || (x >= q.bounds[0] && x < q.bounds[1]) ||
				    (x >= q.bounds[2] && x < q.bounds[3])) slab.push_back(id);
			}
			std::vector<std::vector<CSPrimitives*>> lines(numLines[1]);
			std::vector<bool> haveLine(numLines[1], false);
			std::vector<uint32_t> offsets(1, 0), candidates;
			for (unsigned y = 0; y < numLines[1]; ++y)
			{
				double box[] = {GetDiscLine(0, x ? x-1 : 0), GetDiscLine(0, std::min(x+1,numLines[0]-1)),
				                GetDiscLine(1, y ? y-1 : 0), GetDiscLine(1, std::min(y+1,numLines[1]-1)),
				                GetDiscLine(2, 0), GetDiscLine(2, numLines[2]-1)};
				for (uint32_t id : slab)
				{
					const auto& q = flattened[id];
					bool possible = q.kind == 0;
					// CSXCAD intentionally treats polygon bounds as non-accurate.
					// Bin supported primitives by exact Yee ranges instead, without
					// changing the order or the CPU refinement candidate list.
					for (unsigned n = 0; n < 3 && !possible; ++n)
					{
						unsigned bx = 2*(n == 0), by = 4+2*(n == 1), bz = 8+2*(n == 2);
						possible = x >= q.bounds[bx] && x < q.bounds[bx+1] &&
						           y >= q.bounds[by] && y < q.bounds[by+1] && q.bounds[bz] < q.bounds[bz+1];
					}
					if (possible && primitives[id]->IsInsideBox(box) >= 0) candidates.push_back(id);
				}
				if (candidates.size() > UINT32_MAX) throw std::runtime_error("Metal PEC: too many candidates");
				offsets.push_back(static_cast<uint32_t>(candidates.size()));
			}
			const size_t count = size_t(numLines[1])*numLines[2]*3;
			id<MTLBuffer> rowOffsets = buffer(offsets.data(), offsets.size()*4);
			id<MTLBuffer> rowCandidates = buffer(candidates.data(), candidates.size()*4);
			id<MTLBuffer> output = buffer(nullptr, count*4);
			id<MTLCommandBuffer> command = [queue commandBuffer];
			id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
			[encoder setComputePipelineState:pipeline];
			id<MTLBuffer> buffers[] = {geometry, points, rowOffsets, rowCandidates, coordinates, output};
			for (unsigned i = 0; i < 6; ++i) [encoder setBuffer:buffers[i] offset:0 atIndex:i];
			uint32_t params[] = {numLines[0], numLines[1], numLines[2], x};
			[encoder setBytes:params length:sizeof(params) atIndex:6];
			// Offset the X coordinate without materializing all Yee coordinates.
			[encoder dispatchThreads:MTLSizeMake(numLines[2]*3, numLines[1], 1)
			    threadsPerThreadgroup:MTLSizeMake(pipeline.threadExecutionWidth, 1, 1)];
			[encoder endEncoding]; [command commit]; [command waitUntilCompleted];
			if (command.status == MTLCommandBufferStatusError) throw std::runtime_error("Metal PEC: GPU query failed");
			const uint32_t* winners = static_cast<const uint32_t*>(output.contents);
			for (unsigned y = 0; y < numLines[1]; ++y)
				for (unsigned z = 0; z < numLines[2]; ++z)
					for (unsigned n = 0; n < 3; ++n)
					{
						uint32_t winner = winners[(size_t(y)*numLines[2]+z)*3+n];
						CSPrimitives* prim = winner < primitives.size() ? primitives[winner] : nullptr;
						if (winner == cpuQuery || verify)
						{
							unsigned pos[] = {x,y,z}; double coord[3]; GetYeeCoords(n,pos,coord,false);
							CSPrimitives* reference = nullptr;
							if (!haveLine[y])
							{
								lines[y] = GetPrimitivesBoundBox(x, y, -1, types);
								haveLine[y] = true;
							}
							CSX->GetPropertyByCoordPriority(coord, lines[y], false, &reference);
							if (winner != cpuQuery && prim != reference) throw std::runtime_error("Metal PEC: CPU/GPU winner mismatch");
							prim = reference;
							if (winner == cpuQuery) ++resolved;
						}
						if (prim)
						{
							prim->SetPrimitiveUsed(true);
							if (prim->GetProperty()->GetType() == CSProperties::METAL)
							{
								SetVV(n,x,y,z,0); SetVI(n,x,y,z,0); ++m_Nr_PEC[n];
							}
						}
						++queries;
					}
		}
		CalcPEC_Curves();
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
		std::cout << "Metal PEC: " << queries << " queries, " << resolved << " CPU refinements, "
		          << seconds << " s" << (verify ? " (all winners CPU-verified)" : "") << std::endl;
		return true;
	}
}
