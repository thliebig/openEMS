//! SIMD-accelerated operations for field computations.
//!
//! This module provides vectorized implementations of common FDTD operations
//! using the `wide` crate for portable SIMD with runtime CPU feature detection.

use super::Field3D;
use wide::f32x8;

/// Trait for SIMD-accelerated operations on fields.
pub trait SimdOps {
    /// SIMD-accelerated field update: self += factor * other
    fn simd_fmadd(&mut self, factor: f32, other: &Self);

    /// SIMD-accelerated field update: self = a * self + b * other
    fn simd_axpby(&mut self, a: f32, b: f32, other: &Self);
}

/// Load 8 f32 values from a slice into a SIMD vector.
#[inline]
fn load_f32x8(slice: &[f32]) -> f32x8 {
    let arr: [f32; 8] = slice[..8].try_into().unwrap();
    f32x8::from(arr)
}

/// Store a SIMD vector to a slice.
#[inline]
fn store_f32x8(slice: &mut [f32], v: f32x8) {
    slice[..8].copy_from_slice(v.as_array_ref());
}

impl SimdOps for Field3D {
    fn simd_fmadd(&mut self, factor: f32, other: &Self) {
        let factor_v = f32x8::splat(factor);
        let data = self.as_mut_slice();
        let other_data = other.as_slice();

        // Process 8 elements at a time
        let chunks = data.len() / 8;

        for c in 0..chunks {
            let i = c * 8;
            let a = load_f32x8(&data[i..]);
            let b = load_f32x8(&other_data[i..]);
            let result = a + factor_v * b;
            store_f32x8(&mut data[i..], result);
        }

        // Scalar tail
        for i in (chunks * 8)..data.len() {
            data[i] += factor * other_data[i];
        }
    }

    fn simd_axpby(&mut self, a: f32, b: f32, other: &Self) {
        let a_v = f32x8::splat(a);
        let b_v = f32x8::splat(b);
        let data = self.as_mut_slice();
        let other_data = other.as_slice();

        let chunks = data.len() / 8;

        for c in 0..chunks {
            let i = c * 8;
            let x = load_f32x8(&data[i..]);
            let y = load_f32x8(&other_data[i..]);
            let result = a_v * x + b_v * y;
            store_f32x8(&mut data[i..], result);
        }

        for i in (chunks * 8)..data.len() {
            data[i] = a * data[i] + b * other_data[i];
        }
    }
}

/// Vectorized E-field update for a z-line.
///
/// Computes: E_new = ca * E_old + cb * curl_H
#[inline]
pub fn update_e_line(e_line: &mut [f32], curl_h_line: &[f32], ca: &[f32], cb: &[f32]) {
    debug_assert_eq!(e_line.len(), curl_h_line.len());
    debug_assert_eq!(e_line.len(), ca.len());
    debug_assert_eq!(e_line.len(), cb.len());

    let n = e_line.len();
    let chunks = n / 8;

    // SIMD processing
    for c in 0..chunks {
        let k = c * 8;
        let e = load_f32x8(&e_line[k..]);
        let curl = load_f32x8(&curl_h_line[k..]);
        let ca_v = load_f32x8(&ca[k..]);
        let cb_v = load_f32x8(&cb[k..]);

        let result = ca_v * e + cb_v * curl;
        store_f32x8(&mut e_line[k..], result);
    }

    // Scalar tail
    for k in (chunks * 8)..n {
        e_line[k] = ca[k] * e_line[k] + cb[k] * curl_h_line[k];
    }
}

/// Vectorized H-field update for a z-line.
#[inline]
pub fn update_h_line(h_line: &mut [f32], curl_e_line: &[f32], da: &[f32], db: &[f32]) {
    debug_assert_eq!(h_line.len(), curl_e_line.len());

    let n = h_line.len();
    let chunks = n / 8;

    for c in 0..chunks {
        let k = c * 8;
        let h = load_f32x8(&h_line[k..]);
        let curl = load_f32x8(&curl_e_line[k..]);
        let da_v = load_f32x8(&da[k..]);
        let db_v = load_f32x8(&db[k..]);

        let result = da_v * h + db_v * curl;
        store_f32x8(&mut h_line[k..], result);
    }

    for k in (chunks * 8)..n {
        h_line[k] = da[k] * h_line[k] + db[k] * curl_e_line[k];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrays::Dimensions;

    #[test]
    fn test_simd_fmadd() {
        let dims = Dimensions::new(4, 4, 20);
        let mut a = Field3D::new(dims);
        let b = Field3D::from_vec(dims, vec![1.0; dims.total()]);

        a.fill(2.0);
        a.simd_fmadd(3.0, &b);

        // a = 2.0 + 3.0 * 1.0 = 5.0
        for k in 0..dims.total() {
            assert!((a.as_slice()[k] - 5.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_update_e_line() {
        let n = 100;
        let mut e = vec![1.0f32; n];
        let curl_h = vec![2.0f32; n];
        let ca = vec![0.5f32; n];
        let cb = vec![0.25f32; n];

        update_e_line(&mut e, &curl_h, &ca, &cb);

        // e_new = 0.5 * 1.0 + 0.25 * 2.0 = 0.5 + 0.5 = 1.0
        for v in e.iter() {
            assert!((*v - 1.0).abs() < 1e-6);
        }
    }
}
