// Pipeline constants for grid dimensions
override nx: u32;
override ny: u32;
override nz: u32;

// Fields: [Ex, Ey, Ez] (concatenated)
@group(0) @binding(0) var<storage, read_write> e_field: array<f32>;

// Fields: [Hx, Hy, Hz] (concatenated)
@group(0) @binding(1) var<storage, read_write> h_field: array<f32>;

// E-Coefficients: [Ca_x, Ca_y, Ca_z, Cb_x, Cb_y, Cb_z]
@group(0) @binding(2) var<storage, read> e_coeff: array<f32>;

// H-Coefficients: [Da_x, Da_y, Da_z, Db_x, Db_y, Db_z]
@group(0) @binding(3) var<storage, read> h_coeff: array<f32>;

fn get_idx(x: u32, y: u32, z: u32) -> u32 {
    return x * ny * nz + y * nz + z;
}

// ILP=4: Each thread processes 4 elements along K (Z-axis).
// Workgroup (64, 1, 1) -> processes 256 elements along Z per group.
@compute @workgroup_size(64, 1, 1)
fn update_h(@builtin(global_invocation_id) global_id: vec3<u32>) {
    // k_base is the starting Z-index for this thread
    let k_base = global_id.x * 4u;
    let j = global_id.y;
    let i = global_id.z;

    if (i >= nx || j >= ny || k_base >= nz) {
        return;
    }

    let total = nx * ny * nz;
    let idx_base = get_idx(i, j, k_base);

    // Pre-load data for register reuse along Z-axis (neighbor optimization)
    // We need Ey at k, k+1, k+2, k+3, k+4 for 4 updates.
    // We assume robust buffer access handles out-of-bounds reads with 0 or we clamp.
    // Explicit clamping is safer/faster than relying on robust access overhead.
    
    // Load Ey window
    // Indices for Ey: base, base+1, base+2, base+3, base+4
    // Note: get_idx(i, j, k+1) is just idx + 1 because Z is contiguous
    let ey_idx_base = 1u * total + idx_base;
    
    var ey_vals: array<f32, 5>;
    for (var u = 0u; u < 5u; u++) {
        let k_curr = k_base + u;
        if (k_curr < nz) {
            ey_vals[u] = e_field[ey_idx_base + u];
        } else {
            ey_vals[u] = 0.0;
        }
    }

    // Process 4 elements
    for (var u = 0u; u < 4u; u++) {
        let k = k_base + u;
        if (k >= nz) { break; }
        
        let idx = idx_base + u;

        // Hx update: curl_x = dEz/dy - dEy/dz
        // dez_dy = Ez(i, j+1, k) - Ez(i, j, k)
        // dey_dz = Ey(i, j, k+1) - Ey(i, j, k)
        
        let ez_curr = e_field[2u * total + idx];
        let ez_jp1_val = select(0.0, e_field[2u * total + get_idx(i, j + 1u, k)] - ez_curr, j + 1u < ny);
        let dez_dy = ez_jp1_val;

        // Reuse Ey values from registers
        let dey_dz = ey_vals[u + 1u] - ey_vals[u]; // No need for explicit k+1 check, handled by load loop

        let curl_x = dez_dy - dey_dz;
        let da_x = h_coeff[0u * total + idx];
        let db_x = h_coeff[3u * total + idx];
        let hx_idx = 0u * total + idx;
        h_field[hx_idx] = fma(db_x, curl_x, da_x * h_field[hx_idx]);

        // Hy update: curl_y = dEx/dz - dEz/dx
        // dex_dz = Ex(i, j, k+1) - Ex(i, j, k)
        let ex_curr = e_field[0u * total + idx];
        let ex_kp1 = select(0.0, e_field[0u * total + idx + 1u] - ex_curr, k + 1u < nz);
        let dex_dz = ex_kp1;
        
        let dez_dx = select(0.0, e_field[2u * total + get_idx(i + 1u, j, k)] - ez_curr, i + 1u < nx);
        let curl_y = dex_dz - dez_dx;
        
        let da_y = h_coeff[1u * total + idx];
        let db_y = h_coeff[4u * total + idx];
        let hy_idx = 1u * total + idx;
        h_field[hy_idx] = fma(db_y, curl_y, da_y * h_field[hy_idx]);

        // Hz update: curl_z = dEy/dx - dEx/dy
        let dey_dx = select(0.0, e_field[1u * total + get_idx(i + 1u, j, k)] - ey_vals[u], i + 1u < nx);
        let dex_dy = select(0.0, e_field[0u * total + get_idx(i, j + 1u, k)] - ex_curr, j + 1u < ny);
        let curl_z = dey_dx - dex_dy;
        
        let da_z = h_coeff[2u * total + idx];
        let db_z = h_coeff[5u * total + idx];
        let hz_idx = 2u * total + idx;
        h_field[hz_idx] = fma(db_z, curl_z, da_z * h_field[hz_idx]);
    }
}

@compute @workgroup_size(64, 1, 1)
fn update_e(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let k_base = global_id.x * 4u;
    let j = global_id.y;
    let i = global_id.z;

    if (i < 1u || i >= nx || j < 1u || j >= ny || k_base >= nz) {
        return;
    }

    let total = nx * ny * nz;
    let idx_base = get_idx(i, j, k_base);

    // Pre-load Hy for register reuse
    // We need Hy(i, j, k) and Hy(i, j, k-1)
    // Hy is at index 1.
    // We load Hy[k-1] to Hy[k+3] (5 values)
    // Careful with k=0 (though loop starts at 1, so k-1 is safe if k>=1)
    
    let hy_idx_base = 1u * total + idx_base;
    var hy_vals: array<f32, 5>; // Stores Hy[k-1], Hy[k], Hy[k+1], Hy[k+2], Hy[k+3]
    
    // If k_base=0, k-1 is invalid. But loop condition i<1 || j<1 || k<1 prevents k=0 from processing.
    // So k_base >= 0. But threads 0..3 might have k < 1?
    // Wait, k_base is 0, 4, 8...
    // Thread 0 processes k=0,1,2,3.
    // k=0 is skipped by logic `if (k < 1) continue;` inside loop.
    // So we need valid loads.
    
    for (var u = 0u; u < 5u; u++) {
        // We want Hy at k_base + u - 1
        // indices: -1, 0, 1, 2, 3 relative to k_base
        // Since i,j >= 1, we are safe spatially?
        // k varies.
        // We should check bounds.
        let k_rel = k_base + u; // k_base, k_base+1...
        // But we need k-1.
        
        // Simpler strategy: Load Hy[k]...Hy[k+3] (4 vals) AND Hy[k-1] separately?
        // Or load Hy window from k_base-1 to k_base+3?
        
        // Let's protect the read
        if (k_rel > 0u && k_rel <= nz + 1u) { // Rough check
             let read_k = k_rel - 1u;
             if (read_k < nz) {
                 hy_vals[u] = h_field[hy_idx_base + u - 1u];
             } else {
                 hy_vals[u] = 0.0;
             }
        } else {
             hy_vals[u] = 0.0;
        }
    }
    
    // Process 4 elements
    for (var u = 0u; u < 4u; u++) {
        let k = k_base + u;
        if (k < 1u || k >= nz) { continue; }
        
        let idx = idx_base + u;

        // Ex update: curl_x = (Hz(j) - Hz(j-1)) - (Hy(k) - Hy(k-1))
        let hz_curr = h_field[2u * total + idx];
        let hz_jm1 = h_field[2u * total + get_idx(i, j - 1u, k)];
        
        // Hy(k) is hy_vals[u+1], Hy(k-1) is hy_vals[u]
        let hy_k = hy_vals[u + 1u];
        let hy_km1 = hy_vals[u];
        
        let curl_x = (hz_curr - hz_jm1) - (hy_k - hy_km1);
        
        let ca_x = e_coeff[0u * total + idx];
        let cb_x = e_coeff[3u * total + idx];
        let ex_idx = 0u * total + idx;
        e_field[ex_idx] = fma(cb_x, curl_x, ca_x * e_field[ex_idx]);

        // Ey update: curl_y = (Hx(k) - Hx(k-1)) - (Hz(i) - Hz(i-1))
        let hx_curr = h_field[0u * total + idx];
        let hx_km1 = h_field[0u * total + idx - 1u]; // Safe since k>=1
        let hz_im1 = h_field[2u * total + get_idx(i - 1u, j, k)];
        
        let curl_y = (hx_curr - hx_km1) - (hz_curr - hz_im1);
        
        let ca_y = e_coeff[1u * total + idx];
        let cb_y = e_coeff[4u * total + idx];
        let ey_idx = 1u * total + idx;
        e_field[ey_idx] = fma(cb_y, curl_y, ca_y * e_field[ey_idx]);

        // Ez update: curl_z = (Hy(i) - Hy(i-1)) - (Hx(j) - Hx(j-1))
        let hy_im1 = h_field[1u * total + get_idx(i - 1u, j, k)];
        let hx_jm1 = h_field[0u * total + get_idx(i, j - 1u, k)];
        
        let curl_z = (hy_k - hy_im1) - (hx_curr - hx_jm1);
        
        let ca_z = e_coeff[2u * total + idx];
        let cb_z = e_coeff[5u * total + idx];
        let ez_idx = 2u * total + idx;
        e_field[ez_idx] = fma(cb_z, curl_z, ca_z * e_field[ez_idx]);
    }
}
