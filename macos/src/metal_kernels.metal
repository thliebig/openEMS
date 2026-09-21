//
// openEMS Metal kernels. Compiled at build time into openEMS.metallib and
// embedded in libopenEMS; the runtime never compiles shaders.
//
#include <metal_stdlib>
using namespace metal;
constant bool compressedCoefficients [[function_constant(0)]];
constant bool compressedPML [[function_constant(1)]];
constant bool diamondExcitations [[function_constant(2)]];
constant bool lumpedRLC [[function_constant(3)]];
constant bool diamondUPML [[function_constant(4)]];

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

// One threadgroup owns a space-time diamond for every phase dispatch. The
// host-issued phase boundaries satisfy inter-tile dependencies; barriers here
// order in-place E/source/H/source updates within the tile.
constant uint diamondMaxSteps = 4;

struct DiamondStep
{
	int4 voltage_range; // x begin/end, y begin/end (inclusive)
	int4 current_range;
	uint voltage_source_offset;
	uint voltage_source_count;
	uint current_source_offset;
	uint current_source_count;
	uint rlc_offset;
	uint rlc_count;
	// int4 forces 16-byte alignment, so the record is padded to 64 bytes.
	uint pad0;
	uint pad1;
};

struct DiamondTile
{
	DiamondStep steps[diamondMaxSteps];
};

struct DiamondParams
{
	uint nx;
	uint ny;
	uint nzv;
	uint timestep;
	uint depth;
	float dT_half;
	uint pml_region_count;
	uint pad0;
};

struct DiamondSource
{
	uint field_index;
	float amplitude;
	uint delay;
	uint signal_offset;
	uint signal_length;
	uint period;
};

// One lumped RLC element. The engine packs the trapezoidal state-space
// coefficients so the GPU update is a single read-modify-write of the node
// voltage plus four scalar state words per element.
struct RLCEntry
{
	uint field_index; // scalar float index into the packed voltage buffer
	float dJdV;       // dJ/dV, the implicit self term
	float aV;         // coefficient of Vd[n-1] in the explicit J update
	float aQ;         // coefficient of q[n-1]
	float aJ;         // coefficient of J[n-1]
	float vcd;        // dT/(2*Cd), the node coupling factor
	float vvd;        // 1/(1 + vcd*dJdV)
	float aIl;        // parallel-inductor current coefficient i2v*ilv
};

struct RLCState
{
	float Vd;
	float J;
	float q;
	float Il;
};

// One UPML slab. Coefficients and flux live in the slab's dense local order
// [component][x][y][z]; the diamond kernel applies them per packed lane because
// the four lanes of a float4 span four z positions.
struct PMLRegion
{
	uint sx, sy, sz;
	uint nx, ny, nz;
};

// Argument buffer of per-slab array pointers. Each slab owns six coefficient
// arrays (vv, vvfo, vvfn, ii, iifo, iifn) and two flux arrays (E then H). The
// arrays are the operator's and extension's own allocations, wrapped without a
// copy, so they keep their dense local order.
#define UPML_MAX_REGIONS 6
struct UPMLArgs
{
	const device float* coeff[UPML_MAX_REGIONS * 6];
	device float* flux[UPML_MAX_REGIONS * 2];
};

inline void pml_pre_component(device float4* field, uint fi,
	const device float* self, const device float* old, device float* flux,
	uint off, uint slot, uint nzv, uint sz, uint nz)
{
	const float4 e = field[fi];
	float4 s = 1.0f, o = 0.0f, fv = e;
	for (uint l = 0; l < 4; ++l)
	{
		const uint z = slot + l * nzv;
		if (z >= sz && z < sz + nz)
		{
			const uint zl = z - sz;
			s[l] = self[off + zl];
			o[l] = old[off + zl];
			fv[l] = flux[off + zl];
		}
	}
	const float4 saved = s * e - o * fv;
	field[fi] = fv;
	for (uint l = 0; l < 4; ++l)
	{
		const uint z = slot + l * nzv;
		if (z >= sz && z < sz + nz)
			flux[off + (z - sz)] = saved[l];
	}
}

inline void pml_post_component(device float4* field, uint fi,
	const device float* nw, device float* flux,
	uint off, uint slot, uint nzv, uint sz, uint nz)
{
	const float4 e = field[fi];
	float4 saved = e, nf = 0.0f;
	for (uint l = 0; l < 4; ++l)
	{
		const uint z = slot + l * nzv;
		if (z >= sz && z < sz + nz)
		{
			const uint zl = z - sz;
			saved[l] = flux[off + zl];
			nf[l] = nw[off + zl];
		}
	}
	for (uint l = 0; l < 4; ++l)
	{
		const uint z = slot + l * nzv;
		if (z >= sz && z < sz + nz)
			flux[off + (z - sz)] = e[l];
	}
	field[fi] = saved + nf * e;
}

// Apply one region to all three components of a packed field. Regions are
// disjoint in (x, y), so a cell is touched by at most one slab; callers still
// walk the full list so a nonstandard layout stays correct.
inline void pml_region(device float4* field, uint base, uint x, uint y,
	uint slot, uint nzv, const device PMLRegion* regions, uint r,
	const device UPMLArgs& args, bool voltage, bool post)
{
	const PMLRegion R = regions[r];
	if (x < R.sx || x >= R.sx + R.nx || y < R.sy || y >= R.sy + R.ny)
		return;
	const uint cells = R.nx * R.ny * R.nz;
	const uint cellOff = ((x - R.sx) * R.ny + (y - R.sy)) * R.nz;
	const device float* self = args.coeff[r * 6 + (voltage ? 0 : 3)];
	const device float* old = args.coeff[r * 6 + (voltage ? 1 : 4)];
	const device float* nw = args.coeff[r * 6 + (voltage ? 2 : 5)];
	device float* flux = args.flux[r * 2 + (voltage ? 0 : 1)];
	for (uint n = 0; n < 3; ++n)
	{
		const uint off = n * cells + cellOff;
		if (post)
			pml_post_component(field, base + n, nw, flux, off, slot, nzv, R.sz, R.nz);
		else
			pml_pre_component(field, base + n, self, old, flux, off, slot, nzv, R.sz, R.nz);
	}
}

kernel void update_diamond(
	device float4* volt [[buffer(0)]],
	device float4* curr [[buffer(1)]],
	const device float4* vv [[buffer(2)]],
	const device float4* vi [[buffer(3)]],
	const device float4* ii [[buffer(4)]],
	const device float4* iv [[buffer(5)]],
	constant DiamondParams& p [[buffer(6)]],
	const device DiamondTile* tiles [[buffer(7)]],
	const device ushort* coeffIndex [[buffer(8), function_constant(compressedCoefficients)]],
	const device DiamondSource* sources [[buffer(9), function_constant(diamondExcitations)]],
	const device uint* source_indices [[buffer(10), function_constant(diamondExcitations)]],
	const device float* signal [[buffer(11), function_constant(diamondExcitations)]],
	const device RLCEntry* rlc [[buffer(12), function_constant(lumpedRLC)]],
	device RLCState* rlc_state [[buffer(13), function_constant(lumpedRLC)]],
	const device uint* rlc_indices [[buffer(14), function_constant(lumpedRLC)]],
	const device PMLRegion* pml_regions [[buffer(15), function_constant(diamondUPML)]],
	const device UPMLArgs& pml_args [[buffer(16), function_constant(diamondUPML)]],
	uint tileId [[threadgroup_position_in_grid]],
	uint tid [[thread_index_in_threadgroup]],
	uint threads [[threads_per_threadgroup]])
{
	const uint nzv = p.nzv;
	const uint y_stride = nzv * 3;
	const uint x_stride = p.ny * y_stride;
	// Threads cover whole packed-Z slot groups. Each thread then walks the tile's
	// (x, y) pairs by a fixed stride; the per-step delta is precomputed once, so
	// no runtime integer division appears in the cell loop.
	const uint slots = min(nzv, threads);
	const uint slot0 = tid % slots;
	const uint pair0 = tid / slots;
	const uint groups = threads / slots;
	for (uint timestep = 0; timestep < p.depth; ++timestep)
	{
		const DiamondStep step = tiles[tileId].steps[timestep];
		if (step.voltage_range.x >= 0)
		{
			const uint vx0 = step.voltage_range.x, vy0 = step.voltage_range.z;
			const uint vnx = step.voltage_range.y - vx0 + 1;
			const uint vny = step.voltage_range.w - vy0 + 1;
			const uint pairs = vnx * vny;
			const uint dx = groups / vny, dy = groups % vny;
			for (uint zbase = 0; zbase < nzv; zbase += slots)
			{
				const uint slot = zbase + slot0;
				if (slot >= nzv)
					break;
				uint x = vx0 + pair0 / vny;
				uint y = vy0 + pair0 % vny;
				for (uint pair = pair0; pair < pairs; pair += groups)
				{
					const uint base = x * x_stride + y * y_stride + slot * 3;
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
					if (slot == 0)
					{
						const uint end = base + (nzv - 1) * 3;
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
					if (diamondUPML)
						for (int r = (int)p.pml_region_count - 1; r >= 0; --r)
							pml_region(volt, base, x, y, slot, nzv, pml_regions,
								(uint)r, pml_args, true, false);
					const float4 ex = volt[base];
					const float4 ey = volt[base + 1];
					const float4 ez = volt[base + 2];
					const uint cc = compressedCoefficients ? coeffIndex[base / 3] * 3 : base;
					volt[base] = ex * vv[cc] + vi[cc] * (cz - hz_y - cy + hy_z);
					volt[base + 1] = ey * vv[cc + 1] + vi[cc + 1] * (cx - hx_z - cz + hz_x);
					volt[base + 2] = ez * vv[cc + 2] + vi[cc + 2] * (cy - hy_x - cx + hx_y);
					if (diamondUPML)
						for (uint r = 0; r < p.pml_region_count; ++r)
							pml_region(volt, base, x, y, slot, nzv, pml_regions,
								r, pml_args, true, true);
					y += dy;
					x += dx;
					if (y >= vy0 + vny)
					{
						y -= vny;
						++x;
					}
				}
			}
		}
		threadgroup_barrier(mem_flags::mem_device);
		// Apply2Voltages order on the CPU is by descending extension priority:
		// lumped RLC (default priority) before the excitation (negative), so an
		// edge carrying both sees the RLC node voltage first.
		if ((diamondExcitations || lumpedRLC) && tid == 0)
		{
			if (lumpedRLC)
			{
				device float* field = (device float*)volt;
				for (uint n = step.rlc_offset; n < step.rlc_offset + step.rlc_count; ++n)
				{
					const uint e = rlc_indices[n];
					const RLCEntry entry = rlc[e];
					RLCState s = rlc_state[e];
					const float Vraw = field[entry.field_index];
					s.Il += entry.aIl * s.Vd;
					const float B = entry.aV * s.Vd + entry.aQ * s.q + entry.aJ * s.J;
					const float Vd = entry.vvd * (Vraw - s.Il - entry.vcd * (B + s.J));
					const float J = entry.dJdV * Vd + B;
					s.q += p.dT_half * (J + s.J);
					field[entry.field_index] = Vd;
					s.Vd = Vd;
					s.J = J;
					rlc_state[e] = s;
				}
			}
			if (diamondExcitations)
				for (uint n = step.voltage_source_offset;
				     n < step.voltage_source_offset + step.voltage_source_count; ++n)
				{
					const DiamondSource source = sources[source_indices[n]];
					uint sample = p.timestep + timestep > source.delay ?
						p.timestep + timestep - source.delay : 0;
					sample %= source.period ? source.period : p.timestep + timestep + 1;
					if (sample >= source.signal_length) sample = 0;
					((device float*)volt)[source.field_index] +=
						source.amplitude * signal[source.signal_offset + sample];
				}
		}
		threadgroup_barrier(mem_flags::mem_device);
		if (step.current_range.x >= 0)
		{
			const int stop_x = min(step.current_range.y, (int)p.nx - 2);
			const int stop_y = min(step.current_range.w, (int)p.ny - 2);
			if (stop_x >= step.current_range.x && stop_y >= step.current_range.z)
			{
				const uint cx0 = step.current_range.x, cy0 = step.current_range.z;
				const uint cnx = stop_x - cx0 + 1;
				const uint cny = stop_y - cy0 + 1;
				const uint pairs = cnx * cny;
				const uint dx = groups / cny, dy = groups % cny;
				for (uint zbase = 0; zbase < nzv; zbase += slots)
				{
					const uint slot = zbase + slot0;
					if (slot >= nzv)
						break;
					uint x = cx0 + pair0 / cny;
					uint y = cy0 + pair0 % cny;
					for (uint pair = pair0; pair < pairs; pair += groups)
					{
						const uint base = x * x_stride + y * y_stride + slot * 3;
						const float4 ex = volt[base];
						const float4 ey = volt[base + 1];
						const float4 ez = volt[base + 2];
						const float4 ez_y = volt[base + y_stride + 2];
						const float4 ex_y = volt[base + y_stride];
						const float4 ez_x = volt[base + x_stride + 2];
						const float4 ey_x = volt[base + x_stride + 1];
						float4 ey_z;
						float4 ex_z;
						if (slot + 1 < nzv)
						{
							ey_z = volt[base + 4];
							ex_z = volt[base + 3];
						}
						else
						{
							const uint start = base - slot * 3;
							const float4 ey_start = volt[start + 1];
							const float4 ex_start = volt[start];
							ey_z = float4(ey_start.y, ey_start.z, ey_start.w, 0.0f);
							ex_z = float4(ex_start.y, ex_start.z, ex_start.w, 0.0f);
						}
						if (diamondUPML)
							for (int r = (int)p.pml_region_count - 1; r >= 0; --r)
								pml_region(curr, base, x, y, slot, nzv, pml_regions,
									(uint)r, pml_args, false, false);
						const float4 hx = curr[base];
						const float4 hy = curr[base + 1];
						const float4 hz = curr[base + 2];
						const uint cc = compressedCoefficients ? coeffIndex[base / 3] * 3 : base;
						curr[base] = hx * ii[cc] + iv[cc] * (ez - ez_y - ey + ey_z);
						curr[base + 1] = hy * ii[cc + 1] + iv[cc + 1] * (ex - ex_z - ez + ez_x);
						curr[base + 2] = hz * ii[cc + 2] + iv[cc + 2] * (ey - ey_x - ex + ex_y);
						if (diamondUPML)
							for (uint r = 0; r < p.pml_region_count; ++r)
								pml_region(curr, base, x, y, slot, nzv, pml_regions,
									r, pml_args, false, true);
						y += dy;
						x += dx;
						if (y >= cy0 + cny)
						{
							y -= cny;
							++x;
						}
					}
				}
			}
		}
		threadgroup_barrier(mem_flags::mem_device);
		if (diamondExcitations && tid == 0)
			for (uint n = step.current_source_offset;
			     n < step.current_source_offset + step.current_source_count; ++n)
			{
				const DiamondSource source = sources[source_indices[n]];
				uint sample = p.timestep + timestep > source.delay ?
					p.timestep + timestep - source.delay : 0;
				sample %= source.period ? source.period : p.timestep + timestep + 1;
				if (sample >= source.signal_length) sample = 0;
				((device float*)curr)[source.field_index] +=
					source.amplitude * signal[source.signal_offset + sample];
			}
		threadgroup_barrier(mem_flags::mem_device);
	}
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
