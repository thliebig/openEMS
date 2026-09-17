//
// openEMS Metal kernels. Compiled at build time into openEMS.metallib and
// embedded in libopenEMS; the runtime never compiles shaders.
//
#include <metal_stdlib>
using namespace metal;
constant bool compressedCoefficients [[function_constant(0)]];
constant bool compressedPML [[function_constant(1)]];

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
	const device ushort* coeffIndex [[buffer(5), function_constant(compressedCoefficients)]],
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

struct ExcitationSource
{
	uint fieldIndex;
	float amplitude;
	uint delay;
};

struct ExcitationParams
{
	uint count;
	uint timestep;
	uint signalLength;
	uint period;
};

// Excitations are intentionally applied by one thread. Source lists are sparse,
// and serial application preserves CPU ordering when source entries overlap.
kernel void apply_excitation(
	device float* field [[buffer(0)]],
	const device ExcitationSource* sources [[buffer(1)]],
	const device float* signal [[buffer(2)]],
	constant ExcitationParams& p [[buffer(3)]])
{
	for (uint n = 0; n < p.count; ++n)
	{
		uint sample = p.timestep > sources[n].delay ? p.timestep - sources[n].delay : 0;
		sample %= p.period;
		if (sample >= p.signalLength) sample = 0;
		field[sources[n].fieldIndex] += sources[n].amplitude * signal[sample];
	}
}

kernel void update_currents(
	device float4* curr [[buffer(0)]],
	const device float4* volt [[buffer(1)]],
	const device float4* ii [[buffer(2)]],
	const device float4* iv [[buffer(3)]],
	constant GridParams& p [[buffer(4)]],
	const device ushort* coeffIndex [[buffer(5), function_constant(compressedCoefficients)]],
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

// Plain volt-ADE update for the conducting-sheet model. Each thread owns one
// active Yee edge; the two ADE poles are packed into one float4/float2 record so
// the field is read once for both poles (advance) and the two subtractions run
// in the same order as the CPU reference (apply).
kernel void ade_advance(
	device float* field [[buffer(0)]],
	device float2* state [[buffer(1)]],
	const device float4* coeff [[buffer(2)]],
	const device uint* indices [[buffer(3)]],
	constant uint& count [[buffer(4)]],
	uint gid [[thread_position_in_grid]])
{
	if (gid >= count)
		return;
	const float e = field[indices[gid]];
	float2 s = state[gid];
	const float4 c = coeff[gid];
	s.x = c.x * s.x + c.y * e;
	s.y = c.z * s.y + c.w * e;
	state[gid] = s;
}

kernel void ade_apply(
	device float* field [[buffer(0)]],
	const device float2* state [[buffer(1)]],
	const device uint* indices [[buffer(2)]],
	constant uint& count [[buffer(3)]],
	uint gid [[thread_position_in_grid]])
{
	if (gid >= count)
		return;
	const uint f = indices[gid];
	const float2 s = state[gid];
	float v = field[f];
	v -= s.x;
	v -= s.y;
	field[f] = v;
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

// Auxiliary arrays are reordered once into increasing packed-field addresses.
// Each lane performs only an indexed field access and contiguous auxiliary I/O;
// no coordinate division or component-major passes through the field buffer.
// dispatchThreads supplies exactly the physical component count, including a
// nonuniform final threadgroup; these kernels never address padding elements.
kernel void upml_indexed_pre(
	device float* field [[buffer(0)]],
	device float* flux [[buffer(1)]],
	const device float* self [[buffer(2)]],
	const device float* oldFlux [[buffer(3)]],
	const device uint* indices [[buffer(5)]],
	const device ushort* coeffIndex [[buffer(6), function_constant(compressedPML)]],
	uint gid [[thread_position_in_grid]])
{
	const uint f = indices[gid];
	const uint c = compressedPML ? coeffIndex[gid] : gid;
	const float saved = self[c] * field[f] - oldFlux[c] * flux[gid];
	field[f] = flux[gid];
	flux[gid] = saved;
}

kernel void upml_indexed_post(
	device float* field [[buffer(0)]],
	device float* flux [[buffer(1)]],
	const device float* newFlux [[buffer(2)]],
	const device uint* indices [[buffer(5)]],
	const device ushort* coeffIndex [[buffer(6), function_constant(compressedPML)]],
	uint gid [[thread_position_in_grid]])
{
	const uint f = indices[gid];
	const uint c = compressedPML ? coeffIndex[gid] : gid;
	const float saved = flux[gid];
	flux[gid] = field[f];
	field[f] = saved + newFlux[c] * flux[gid];
}

// PEC geometry pass (host test and shader share metal_predicates.h).
#include "metal_predicates.h"

struct Primitive { uint bounds[12]; uint kind, normal, first, count; };
struct Vertex { float xh, xl, yh, yl; };
struct Cylinder { float p0[3], radius, p1[3], shell; };

kernel void pec_mask(const device Primitive* prims [[buffer(0)]],
                     const device Vertex* vertices [[buffer(1)]],
                     const device uint* offsets [[buffer(2)]],
                     const device uint* candidates [[buffer(3)]],
                     const device float* grid [[buffer(4)]],
                     device uint* winners [[buffer(5)]],
                     constant uint4& p [[buffer(6)]],
                     const device Cylinder* cyls [[buffer(7)]],
                     const device float2* gridD [[buffer(8)]],
                     constant uint& dualMesh [[buffer(9)]],
                     uint3 gid [[thread_position_in_grid]])
{
	if (gid.x >= p.z * 3 || gid.y >= p.y) return;
	uint n = gid.x % 3, z = gid.x / 3, y = gid.y, x = p.w;
	uint line = y;
	uint output = line * p.z * 3 + gid.x;
	// dualMesh mirrors Operator::GetYeeCoords(n,pos,coord,dualMesh): the component
	// axis takes the dual offset for the primal/voltage query and the primal
	// offset for the dual/current query, the two transverse axes the opposite.
	uint dn = dualMesh ? 0u : 1u;
	uint dt = dualMesh ? 1u : 0u;
	uint g0 = 2*x + (n == 0 ? dn : dt);
	uint g1 = 2*p.x + 2*y + (n == 1 ? dn : dt);
	uint g2 = 2*(p.x+p.y) + 2*z + (n == 2 ? dn : dt);
	float c[3] = {grid[g0], grid[g1], grid[g2]};
	uint gidx[3] = {g0, g1, g2};
	mp_df cd[3] = {{gridD[gidx[0]].x, gridD[gidx[0]].y},
	               {gridD[gidx[1]].x, gridD[gidx[1]].y},
	               {gridD[gidx[2]].x, gridD[gidx[2]].y}};
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
			uint b = a*4 + 2*((a == n) ? dn : dt);
			outside |= pos[a] < q.bounds[b] || pos[a] >= q.bounds[b+1];
		}
		if (outside) continue;
		if (q.kind == 1) { result = id; break; }
		if (q.kind == 3 || q.kind == 4)
		{
			Cylinder cy = cyls[q.first];
			float3 a = float3(cy.p0[0], cy.p0[1], cy.p0[2]);
			float3 b = float3(cy.p1[0], cy.p1[1], cy.p1[2]);
			float3 pp = float3(c[0], c[1], c[2]);
			float3 ab = b - a;
			float ab2 = dot(ab, ab);
			float scale = max(1.0f, max(abs(cy.radius), abs(cy.shell)));
			scale = max(scale, max(max(abs(pp.x), abs(pp.y)), abs(pp.z)));
			scale = max(scale, max(max(abs(a.x), abs(a.y)), max(abs(b.x), max(abs(b.y), abs(b.z)))));
			float e = 256.0f * FLT_EPSILON * scale;
			float te = 256.0f * FLT_EPSILON;
			if (ab2 <= e*e)
			{
				uncertain = true;
			}
			else
			{
				float t = dot(pp - a, ab) / ab2;
				if (t < -te || t > 1.0f + te)
				{
					// outside the axis segment, not this primitive
				}
				else if (t < te || t > 1.0f - te)
				{
					uncertain = true;
				}
				else
				{
					float d = length(pp - (a + t*ab));
					if (q.kind == 3)
					{
						if (d <= cy.radius - e) { result = id; break; }
						if (d <= cy.radius + e) uncertain = true;
					}
					else
					{
						float lower = cy.radius - 0.5f*cy.shell;
						float upper = cy.radius + 0.5f*cy.shell;
						if (d < lower - e || d > upper + e) { }
						else if (d <= lower + e || d >= upper - e) uncertain = true;
						else { result = id; break; }
					}
				}
			}
			if (uncertain) { result = 0xfffffffeu; break; }
			continue;
		}
		// Mirror CSPrimPolygon::IsInside with exact (double-float) predicates.
		mp_df px = cd[(q.normal+1)%3];
		mp_df py = cd[(q.normal+2)%3];
		uint vlast = q.first + q.count - 1;
		mp_df x1 = {vertices[vlast].xh, vertices[vlast].xl};
		mp_df y1 = {vertices[vlast].yh, vertices[vlast].yl};
		int winding = 0;
		bool onedge = false;
		int startover = mp_df_ge(y1, py);
		for (uint j = 0; j < q.count; ++j)
		{
			uint vj = q.first + j;
			mp_df x2 = {vertices[vj].xh, vertices[vj].xl};
			mp_df y2 = {vertices[vj].yh, vertices[vj].yl};
			// Exact axis-aligned on-edge tests.
			if (mp_df_eq(x2, x1) && mp_df_eq(x1, px) &&
			    ((!mp_df_ge(py, y1) && !mp_df_ge(y2, py)) ||
			     (!mp_df_ge(y1, py) && !mp_df_ge(py, y2))))
			{ onedge = true; break; }
			if (mp_df_eq(y2, y1) && mp_df_eq(y1, py) &&
			    ((!mp_df_ge(px, x1) && !mp_df_ge(x2, px)) ||
			     (!mp_df_ge(x1, px) && !mp_df_ge(px, x2))))
			{ onedge = true; break; }
			int endover = mp_df_ge(y2, py);
			if (startover != endover)
			{
				int s = mp_orient2d_sign(x1, y1, x2, y2, px, py);
				if (s == 0) { uncertain = true; break; }
				// CSXCAD: (y2-py)*(x2-x1) <= (y2-y1)*(x2-px), i.e. orient2d >= 0.
				if (s > 0) { if (endover) ++winding; }
				else { if (!endover) --winding; }
			}
			startover = endover;
			x1 = x2; y1 = y2;
		}
		if (uncertain) { result = 0xfffffffeu; break; }
		if (onedge || winding != 0) { result = id; break; }
		continue;
	}
	winners[output] = result;
}
