/*
 * Copyright (C) 2026 openEMS contributors
 * SPDX-License-Identifier: GPL-3.0-or-later
 */
#include "operator_metal.h"
#include "ContinuousStructure.h"
#include "CSPrimPolygon.h"
#include "CSPrimLinPoly.h"
#include "CSPrimCylinder.h"
#include "CSPrimCylindricalShell.h"
#include "CSPropConductingSheet.h"
#include "extensions/operator_ext_conductingsheet.h"
#include "extensions/operator_ext_lorentzmaterial.h"
#include "metal_library.h"

#import <Foundation/Foundation.h>

#include <string>
#import <Metal/Metal.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <typeinfo>
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
struct PecVertex { float xh, xl, yh, yl; };
struct PecCoord { float hi, lo; };
struct PecCylinder { float p0[3], radius, p1[3], shell; };
const uint32_t cpuQuery = UINT32_MAX - 1;

// Class of EC-consuming extension that owns a resolved winner: 0 = none,
// 1 = conducting sheet, 2 = dispersive (Lorentz/Debye).
uint8_t ClassifyGeometryWinner(CSPrimitives* prim)
{
	if (!prim)
		return 0;
	CSProperties* prop = prim->GetProperty();
	if (!prop)
		return 0;
	if (dynamic_cast<CSPropConductingSheet*>(prop))
		return 1;
	if (prop->ToLorentzMaterial() || prop->ToDebyeMaterial())
		return 2;
	return 0;
}
}

bool MetalDeviceAvailable(std::string& reason)
{
	@autoreleasepool
	{
		id<MTLDevice> device = MTLCreateSystemDefaultDevice();
		if (device)
			return true;
	}
	reason = "no Metal GPU device available";
	return false;
}

bool Operator_Metal::CalcPEC()
{
	const auto started = std::chrono::steady_clock::now();
	const char* setting = std::getenv("OPENEMS_METAL_PEC");
	m_geoConductingSheet.clear();
	m_geoDispersivePrimal.clear();
	m_geoDispersiveDual.clear();
	m_geoWinnersValid = false;
	auto cpuFallback = [&](const char* why) {
		std::cerr << "Metal PEC: " << why << "; using CPU PEC mapping" << std::endl;
		return Operator::CalcPEC();
	};
	// --engine=metal enables PEC mapping too; retain a diagnostic CPU override.
	if (setting && setting[0] == '0')
	{
		const bool result = cpuFallback("selected by OPENEMS_METAL_PEC=0");
		std::cout << "CPU PEC: " << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count() << " s" << std::endl;
		return result;
	}
	if (m_MeshType != CARTESIAN)
		return cpuFallback("non-Cartesian mesh is not supported");
	const bool verify = setting && std::string(setting) == "verify";
	@autoreleasepool
	{
		id<MTLDevice> device = MTLCreateSystemDefaultDevice();
		id<MTLCommandQueue> queue = [device newCommandQueue];
		NSError* error = nil;
		id<MTLLibrary> library = OpenEMSMetalLibrary(device, &error);
		id<MTLFunction> function = library ? [library newFunctionWithName:@"pec_mask"] : nil;
		id<MTLComputePipelineState> pipeline = function ? [device newComputePipelineStateWithFunction:function error:&error] : nil;
		if (!queue || !pipeline)
		{
			std::cerr << "Metal PEC: unavailable; using CPU: " << (error ? error.localizedDescription.UTF8String : "no device") << std::endl;
			return Operator::CalcPEC();
		}
		auto buffer = [device](const void* data, size_t bytes) {
			// Metal rejects zero-length buffers; all geometry element types are at
			// least one float wide, so a float-sized dummy is a safe minimum.
			id<MTLBuffer> b = data && bytes ? [device newBufferWithBytes:data length:bytes options:MTLResourceStorageModeShared]
			    : [device newBufferWithLength:std::max<size_t>(bytes, sizeof(float)) options:MTLResourceStorageModeShared];
			if (!b) throw std::runtime_error("Metal PEC: buffer allocation failed");
			return b;
		};
		const auto types = static_cast<CSProperties::PropertyType>(CSProperties::MATERIAL | CSProperties::METAL);
		const auto primitives = CSX->GetAllPrimitives(true, types);
		if (primitives.size() >= cpuQuery) return cpuFallback("primitive count exceeds the GPU index range");
		// Convert bounds to exact Yee-index ranges on the CPU. In particular,
		// zero-thickness sheets must not become a tolerance-thickened volume.
		std::vector<double> axes[3][2];
		std::vector<float> grid;
		std::vector<PecCoord> gridD;
		for (int a = 0; a < 3; ++a)
			for (unsigned i = 0; i < numLines[a]; ++i)
				for (int dual = 0; dual < 2; ++dual)
				{
					double coord = GetDiscLine(a, i, dual);
					if (!std::isfinite(coord) || std::abs(coord) > 1e10) return cpuFallback("non-finite or out-of-range grid coordinates");
					axes[a][dual].push_back(coord);
					grid.push_back(static_cast<float>(coord));
					PecCoord gc;
					gc.hi = static_cast<float>(coord);
					gc.lo = static_cast<float>(coord - static_cast<double>(gc.hi));
					gridD.push_back(gc);
				}
		for (int a = 0; a < 3; ++a)
			for (int dual = 0; dual < 2; ++dual)
				if (!std::is_sorted(axes[a][dual].begin(), axes[a][dual].end())) return cpuFallback("grid lines are not sorted");
		std::vector<PecPrimitive> flattened;
		std::vector<PecVertex> vertices;
		std::vector<PecCylinder> cylinders;
		size_t unsupported = 0;
		for (auto* prim : primitives)
		{
			PecPrimitive q{};
			const int type = prim->GetType();
			const bool cartesian_prim =
			    !prim->HasTransform() && prim->GetCoordInputType() == CARTESIAN &&
			    (prim->GetCoordinateSystem() == CARTESIAN || prim->GetCoordinateSystem() == UNDEFINED_CS);
			const bool flat_supported = cartesian_prim &&
			    (type == CSPrimitives::BOX || type == CSPrimitives::POLYGON || type == CSPrimitives::LINPOLY);
			const bool cylinder_supported = cartesian_prim &&
			    (type == CSPrimitives::CYLINDER || type == CSPrimitives::CYLINDRICALSHELL);
			if (flat_supported || cylinder_supported)
			{
				double bounds[6];
				prim->GetBoundBox(bounds);
				if (flat_supported)
					q.kind = type == CSPrimitives::BOX ? 1 : 2;
				else
					q.kind = type == CSPrimitives::CYLINDER ? 3 : 4;
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
						return cpuFallback("polygon has too many vertices");
					q.normal = polygon->GetNormDir();
					q.first = static_cast<uint32_t>(vertices.size());
					q.count = static_cast<uint32_t>(polygon->GetQtyCoords());
					if (!q.count || q.normal > 2) q.kind = 0;
					for (uint32_t j = 0; j < q.count; ++j)
					{
						double vx = polygon->GetCoord(2*j), vy = polygon->GetCoord(2*j+1);
						if (!std::isfinite(vx) || !std::isfinite(vy) || std::abs(vx) > 1e10 || std::abs(vy) > 1e10) q.kind = 0;
						PecVertex vert;
						vert.xh = static_cast<float>(vx);
						vert.xl = static_cast<float>(vx - static_cast<double>(vert.xh));
						vert.yh = static_cast<float>(vy);
						vert.yl = static_cast<float>(vy - static_cast<double>(vert.yh));
						vertices.push_back(vert);
					}
				}
				else if (q.kind == 3 || q.kind == 4)
				{
					// Cylinders and cylindrical shells share the axis/radius test; a
					// shell additionally bounds the distance to its radius. Vias are
					// commonly cylinders, so both are mapped on the GPU.
					auto* cylinder = static_cast<CSPrimCylinder*>(prim);
					const double* start = cylinder->GetAxisStartCoord()->GetCartesianCoords();
					const double* stop = cylinder->GetAxisStopCoord()->GetCartesianCoords();
					const double radius = cylinder->GetRadius();
					const double shell = q.kind == 4
					    ? static_cast<CSPrimCylindricalShell*>(prim)->GetShellWidth() : 0.0;
					bool valid = std::isfinite(radius) && std::isfinite(shell) && radius >= 0 && shell >= 0;
					for (int i = 0; i < 3; ++i)
						valid = valid && std::isfinite(start[i]) && std::isfinite(stop[i]) &&
						        std::abs(start[i]) <= 1e10 && std::abs(stop[i]) <= 1e10;
					const double dx = stop[0]-start[0], dy = stop[1]-start[1], dz = stop[2]-start[2];
					valid = valid && (dx*dx + dy*dy + dz*dz) > 0.0;
					if (!valid) q.kind = 0;
					else
					{
						PecCylinder c{};
						c.p0[0] = static_cast<float>(start[0]);
						c.p0[1] = static_cast<float>(start[1]);
						c.p0[2] = static_cast<float>(start[2]);
						c.radius = static_cast<float>(radius);
						c.p1[0] = static_cast<float>(stop[0]);
						c.p1[1] = static_cast<float>(stop[1]);
						c.p1[2] = static_cast<float>(stop[2]);
						c.shell = static_cast<float>(shell);
						q.first = static_cast<uint32_t>(cylinders.size());
						cylinders.push_back(c);
					}
				}
			}
			if (!q.kind) ++unsupported;
			flattened.push_back(q);
		}
		if (unsupported)
			std::cout << "Metal PEC: " << unsupported << " unsupported primitives; affected queries use CPU" << std::endl;
		// EC-consuming extensions (lossy conductor / dispersive dielectric) need the
		// same MATERIAL|METAL winner at each Yee component. Record theirs during this
		// pass so they do not re-collect and re-sort every primitive per (x,y) row.
		bool needConductingSheet = false, needDispersive = false;
		for (auto* ext : m_Op_exts)
		{
			if (typeid(*ext) == typeid(Operator_Ext_ConductingSheet)) needConductingSheet = true;
			if (typeid(*ext) == typeid(Operator_Ext_LorentzMaterial)) needDispersive = true;
		}
		const bool recordWinners = needConductingSheet || needDispersive;
		// Classify primitives once so the per-winner test is an array read, not RTTI.
		std::vector<uint8_t> primClass(primitives.size(), 0);
		for (size_t id = 0; id < primitives.size(); ++id)
			primClass[id] = ClassifyGeometryWinner(primitives[id]);
		id<MTLBuffer> geometry = buffer(flattened.data(), flattened.size()*sizeof(PecPrimitive));
		id<MTLBuffer> points = buffer(vertices.data(), vertices.size()*sizeof(PecVertex));
		id<MTLBuffer> cylinderGeometry = buffer(cylinders.data(), cylinders.size()*sizeof(PecCylinder));
		id<MTLBuffer> coordinates = buffer(grid.data(), grid.size()*sizeof(float));
		id<MTLBuffer> coordinatesD = buffer(gridD.data(), gridD.size()*sizeof(PecCoord));
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
			std::vector<std::vector<CSPrimitives*>> refLines(numLines[1]);
			std::vector<bool> haveRef(numLines[1], false);
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
					if (!possible)
					{
						// A component query may use the primal or the dual line set on
						// each axis (primal and dual dispatches share this list), so
						// include a primitive that could contain the cell on either.
						// The shader still decides against the exact bounds.
						const bool in_x = (x >= q.bounds[0] && x < q.bounds[1]) ||
						                  (x >= q.bounds[2] && x < q.bounds[3]);
						const bool in_y = (y >= q.bounds[4] && y < q.bounds[5]) ||
						                  (y >= q.bounds[6] && y < q.bounds[7]);
						const bool z_any = q.bounds[8] < q.bounds[9] || q.bounds[10] < q.bounds[11];
						possible = in_x && in_y && z_any;
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
			// Buffer 6 is the inline params struct set below, so the cylinder and
			// high-precision grid buffers start at index 7.
			id<MTLBuffer> buffers[] = {geometry, points, rowOffsets, rowCandidates, coordinates, output};
			for (unsigned i = 0; i < 6; ++i) [encoder setBuffer:buffers[i] offset:0 atIndex:i];
			[encoder setBuffer:cylinderGeometry offset:0 atIndex:7];
			[encoder setBuffer:coordinatesD offset:0 atIndex:8];
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
								// The per-row candidate list already is the bounding-box filtered,
								// CSXCAD priority-ordered primitive list, so reuse it instead of
								// re-fetching and re-sorting every primitive for this row.
								lines[y].clear();
								lines[y].reserve(offsets[y+1]-offsets[y]);
								for (uint32_t k = offsets[y]; k < offsets[y+1]; ++k)
									lines[y].push_back(primitives[candidates[k]]);
								haveLine[y] = true;
							}
							CSX->GetPropertyByCoordPriority(coord, lines[y], false, &reference);
							if (verify)
							{
								// Cross-check against the full bounding-box list. The row
								// prefilter may only drop primitives that cannot contain this
								// coordinate, so the resolved winner must be identical.
								if (!haveRef[y])
								{
									double box[] = {GetDiscLine(0, x ? x-1 : 0), GetDiscLine(0, std::min(x+1,numLines[0]-1)),
									                GetDiscLine(1, y ? y-1 : 0), GetDiscLine(1, std::min(y+1,numLines[1]-1)),
									                GetDiscLine(2, 0), GetDiscLine(2, numLines[2]-1)};
									for (auto* candidate : primitives)
										if (candidate->IsInsideBox(box) >= 0) refLines[y].push_back(candidate);
									haveRef[y] = true;
								}
								CSPrimitives* authority = nullptr;
								CSX->GetPropertyByCoordPriority(coord, refLines[y], false, &authority);
								if (reference != authority || (winner != cpuQuery && prim != authority))
									throw std::runtime_error("Metal PEC: CPU/GPU winner mismatch");
							}
							prim = reference;
							if (winner == cpuQuery) ++resolved;
						}
						if (prim)
						{
							if (recordWinners)
							{
								uint8_t cls = 0;
								if (!verify && winner < primitives.size())
									cls = primClass[winner];
								else
									cls = ClassifyGeometryWinner(prim);
								if (cls == 1 && needConductingSheet)
									m_geoConductingSheet.push_back({x,y,z,static_cast<unsigned char>(n),prim});
								else if (cls == 2 && needDispersive)
									m_geoDispersivePrimal.push_back({x,y,z,static_cast<unsigned char>(n),prim});
							}
							prim->SetPrimitiveUsed(true);
							if (prim->GetProperty()->GetType() == CSProperties::METAL)
							{
								SetVV(n,x,y,z,0); SetVI(n,x,y,z,0); ++m_Nr_PEC[n];
							}
						}
						++queries;
					}
			// Dispersive current (magnetic) terms query the dual grid. Resolve a second
			// winner set for the same slab, sharing candidates, only when needed.
			if (needDispersive)
			{
				id<MTLBuffer> outputD = buffer(nullptr, count*4);
				id<MTLCommandBuffer> commandD = [queue commandBuffer];
				id<MTLComputeCommandEncoder> encoderD = [commandD computeCommandEncoder];
				[encoderD setComputePipelineState:pipeline];
				id<MTLBuffer> buffersD[] = {geometry, points, rowOffsets, rowCandidates, coordinates, outputD};
				for (unsigned i = 0; i < 6; ++i) [encoderD setBuffer:buffersD[i] offset:0 atIndex:i];
				[encoderD setBuffer:cylinderGeometry offset:0 atIndex:7];
				[encoderD setBuffer:coordinatesD offset:0 atIndex:8];
				[encoderD setBytes:params length:sizeof(params) atIndex:6];
				uint32_t dualFlag = 1;
				[encoderD setBytes:&dualFlag length:sizeof(dualFlag) atIndex:9];
				[encoderD dispatchThreads:MTLSizeMake(numLines[2]*3, numLines[1], 1)
				    threadsPerThreadgroup:MTLSizeMake(pipeline.threadExecutionWidth, 1, 1)];
				[encoderD endEncoding]; [commandD commit]; [commandD waitUntilCompleted];
				if (commandD.status == MTLCommandBufferStatusError) throw std::runtime_error("Metal PEC: GPU dual query failed");
				const uint32_t* winnersD = static_cast<const uint32_t*>(outputD.contents);
				for (unsigned y = 0; y < numLines[1]; ++y)
					for (unsigned z = 0; z < numLines[2]; ++z)
						for (unsigned n = 0; n < 3; ++n)
						{
							uint32_t winner = winnersD[(size_t(y)*numLines[2]+z)*3+n];
							CSPrimitives* prim = winner < primitives.size() ? primitives[winner] : nullptr;
							uint8_t cls = winner < primitives.size() ? primClass[winner] : 0;
							if (winner == cpuQuery)
							{
								unsigned pos[] = {x,y,z}; double coord[3];
								if (GetYeeCoords(n,pos,coord,true)==false) continue;
								if (!haveLine[y])
								{
									lines[y].clear();
									lines[y].reserve(offsets[y+1]-offsets[y]);
									for (uint32_t k = offsets[y]; k < offsets[y+1]; ++k)
										lines[y].push_back(primitives[candidates[k]]);
									haveLine[y] = true;
								}
								CSX->GetPropertyByCoordPriority(coord, lines[y], false, &prim);
								cls = ClassifyGeometryWinner(prim);
							}
							if (cls == 2 && prim)
								m_geoDispersiveDual.push_back({x,y,z,static_cast<unsigned char>(n),prim});
						}
			}
		}
		m_geoWinnersValid = true;
		CalcPEC_Curves();
		const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
		std::cout << "Metal PEC: " << queries << " queries, " << resolved << " CPU refinements, "
		          << seconds << " s" << (verify ? " (all winners CPU-verified)" : "") << std::endl;
		return true;
	}
}
