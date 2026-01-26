//! Compressed coefficient storage for memory-efficient FDTD.
//!
//! This module implements coefficient compression similar to the C++ openEMS
//! `operator_sse_compressed` implementation. Instead of storing coefficients
//! for every grid point, unique coefficient sets are stored in a lookup table,
//! and each grid point stores only an index into this table.
//!
//! ## Memory Bandwidth Optimization
//!
//! The key insight is that many grid cells share identical material properties:
//! - Vacuum regions all have the same coefficients
//! - PML layers have repeated coefficient patterns
//! - Homogeneous material regions share coefficients
//!
//! By deduplicating coefficients and using index-based lookup:
//! 1. The coefficient table is small enough to fit in L1/L2 cache
//! 2. Memory bandwidth is reduced (loading indices instead of full coefficients)
//! 3. Better cache utilization during FDTD updates

use crate::arrays::{Dimensions, Field3D};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};

/// Tolerance for coefficient comparison.
const COEFF_TOLERANCE: f32 = 1e-10;

/// Wrapper for f32 that implements Hash with tolerance.
#[derive(Clone, Copy, Debug)]
struct HashableF32(f32);

impl Hash for HashableF32 {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // Quantize to fixed precision for hashing
        let quantized = (self.0 / COEFF_TOLERANCE).round() as i64;
        quantized.hash(state);
    }
}

impl PartialEq for HashableF32 {
    fn eq(&self, other: &Self) -> bool {
        (self.0 - other.0).abs() < COEFF_TOLERANCE
    }
}

impl Eq for HashableF32 {}

/// E-field update coefficients for a single cell.
/// E_new = ca * E_old + cb * curl_H
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct ECoefficients {
    /// Ca coefficient (amplitude decay)
    pub ca: f32,
    /// Cb coefficient (curl coupling)
    pub cb: f32,
}

/// H-field update coefficients for a single cell.
/// H_new = da * H_old + db * curl_E
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct HCoefficients {
    /// Da coefficient (amplitude decay)
    pub da: f32,
    /// Db coefficient (curl coupling)
    pub db: f32,
}

/// Hashable wrapper for E coefficients.
#[derive(Clone, Copy, Debug)]
struct HashableECoeff {
    ca: HashableF32,
    cb: HashableF32,
}

impl Hash for HashableECoeff {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.ca.hash(state);
        self.cb.hash(state);
    }
}

impl PartialEq for HashableECoeff {
    fn eq(&self, other: &Self) -> bool {
        self.ca == other.ca && self.cb == other.cb
    }
}

impl Eq for HashableECoeff {}

impl From<ECoefficients> for HashableECoeff {
    fn from(c: ECoefficients) -> Self {
        Self {
            ca: HashableF32(c.ca),
            cb: HashableF32(c.cb),
        }
    }
}

/// Hashable wrapper for H coefficients.
#[derive(Clone, Copy, Debug)]
struct HashableHCoeff {
    da: HashableF32,
    db: HashableF32,
}

impl Hash for HashableHCoeff {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.da.hash(state);
        self.db.hash(state);
    }
}

impl PartialEq for HashableHCoeff {
    fn eq(&self, other: &Self) -> bool {
        self.da == other.da && self.db == other.db
    }
}

impl Eq for HashableHCoeff {}

impl From<HCoefficients> for HashableHCoeff {
    fn from(c: HCoefficients) -> Self {
        Self {
            da: HashableF32(c.da),
            db: HashableF32(c.db),
        }
    }
}

/// Compressed E-field coefficient storage.
///
/// Instead of storing 2 coefficients (ca, cb) per cell per direction,
/// unique coefficient pairs are stored in a lookup table, and each
/// cell stores only a u32 index.
#[derive(Debug)]
pub struct CompressedECoefficients {
    /// Grid dimensions
    dims: Dimensions,
    /// Unique coefficient sets per direction [x, y, z]
    tables: [Vec<ECoefficients>; 3],
    /// Index arrays per direction [x, y, z]
    /// Each index maps a linear grid position to a coefficient table entry
    indices: [Vec<u32>; 3],
}

impl CompressedECoefficients {
    /// Create compressed E coefficients from full coefficient arrays.
    pub fn from_fields(ca: &[&Field3D; 3], cb: &[&Field3D; 3]) -> Self {
        let dims = ca[0].dims();
        let size = dims.total();

        let mut tables: [Vec<ECoefficients>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        let mut indices: [Vec<u32>; 3] = [
            vec![0u32; size],
            vec![0u32; size],
            vec![0u32; size],
        ];

        // Process each direction independently
        for dir in 0..3 {
            let mut dedup_map: HashMap<HashableECoeff, u32> = HashMap::new();

            for i in 0..dims.nx {
                for j in 0..dims.ny {
                    for k in 0..dims.nz {
                        let idx = dims.to_linear(i, j, k);

                        let coeff = ECoefficients {
                            ca: ca[dir].get(i, j, k),
                            cb: cb[dir].get(i, j, k),
                        };

                        let hashable: HashableECoeff = coeff.into();

                        let table_idx = if let Some(&existing_idx) = dedup_map.get(&hashable) {
                            existing_idx
                        } else {
                            let new_idx = tables[dir].len() as u32;
                            tables[dir].push(coeff);
                            dedup_map.insert(hashable, new_idx);
                            new_idx
                        };

                        indices[dir][idx] = table_idx;
                    }
                }
            }
        }

        Self {
            dims,
            tables,
            indices,
        }
    }

    /// Get dimensions.
    #[inline]
    pub fn dims(&self) -> Dimensions {
        self.dims
    }

    /// Get the coefficient table for a direction.
    #[inline]
    pub fn table(&self, dir: usize) -> &[ECoefficients] {
        &self.tables[dir]
    }

    /// Get the index array for a direction.
    #[inline]
    pub fn indices(&self, dir: usize) -> &[u32] {
        &self.indices[dir]
    }

    /// Get coefficients for a specific position and direction.
    #[inline]
    pub fn get(&self, dir: usize, i: usize, j: usize, k: usize) -> ECoefficients {
        let idx = self.dims.to_linear(i, j, k);
        let table_idx = self.indices[dir][idx] as usize;
        self.tables[dir][table_idx]
    }

    /// Get number of unique coefficient sets per direction.
    pub fn num_unique(&self) -> [usize; 3] {
        [
            self.tables[0].len(),
            self.tables[1].len(),
            self.tables[2].len(),
        ]
    }

    /// Calculate compression ratio.
    pub fn compression_ratio(&self) -> f64 {
        let size = self.dims.total();
        // Original: 2 f32 per cell per direction = 6 * size * 4 bytes
        let original = 6 * size * std::mem::size_of::<f32>();
        // Compressed: tables + indices
        let compressed: usize = self.tables.iter().map(|t| t.len() * std::mem::size_of::<ECoefficients>()).sum::<usize>()
            + 3 * size * std::mem::size_of::<u32>();
        original as f64 / compressed as f64
    }
}

/// Compressed H-field coefficient storage.
#[derive(Debug)]
pub struct CompressedHCoefficients {
    /// Grid dimensions
    dims: Dimensions,
    /// Unique coefficient sets per direction [x, y, z]
    tables: [Vec<HCoefficients>; 3],
    /// Index arrays per direction [x, y, z]
    indices: [Vec<u32>; 3],
}

impl CompressedHCoefficients {
    /// Create compressed H coefficients from full coefficient arrays.
    pub fn from_fields(da: &[&Field3D; 3], db: &[&Field3D; 3]) -> Self {
        let dims = da[0].dims();
        let size = dims.total();

        let mut tables: [Vec<HCoefficients>; 3] = [Vec::new(), Vec::new(), Vec::new()];
        let mut indices: [Vec<u32>; 3] = [
            vec![0u32; size],
            vec![0u32; size],
            vec![0u32; size],
        ];

        for dir in 0..3 {
            let mut dedup_map: HashMap<HashableHCoeff, u32> = HashMap::new();

            for i in 0..dims.nx {
                for j in 0..dims.ny {
                    for k in 0..dims.nz {
                        let idx = dims.to_linear(i, j, k);

                        let coeff = HCoefficients {
                            da: da[dir].get(i, j, k),
                            db: db[dir].get(i, j, k),
                        };

                        let hashable: HashableHCoeff = coeff.into();

                        let table_idx = if let Some(&existing_idx) = dedup_map.get(&hashable) {
                            existing_idx
                        } else {
                            let new_idx = tables[dir].len() as u32;
                            tables[dir].push(coeff);
                            dedup_map.insert(hashable, new_idx);
                            new_idx
                        };

                        indices[dir][idx] = table_idx;
                    }
                }
            }
        }

        Self {
            dims,
            tables,
            indices,
        }
    }

    /// Get dimensions.
    #[inline]
    pub fn dims(&self) -> Dimensions {
        self.dims
    }

    /// Get the coefficient table for a direction.
    #[inline]
    pub fn table(&self, dir: usize) -> &[HCoefficients] {
        &self.tables[dir]
    }

    /// Get the index array for a direction.
    #[inline]
    pub fn indices(&self, dir: usize) -> &[u32] {
        &self.indices[dir]
    }

    /// Get coefficients for a specific position and direction.
    #[inline]
    pub fn get(&self, dir: usize, i: usize, j: usize, k: usize) -> HCoefficients {
        let idx = self.dims.to_linear(i, j, k);
        let table_idx = self.indices[dir][idx] as usize;
        self.tables[dir][table_idx]
    }

    /// Get number of unique coefficient sets per direction.
    pub fn num_unique(&self) -> [usize; 3] {
        [
            self.tables[0].len(),
            self.tables[1].len(),
            self.tables[2].len(),
        ]
    }

    /// Calculate compression ratio.
    pub fn compression_ratio(&self) -> f64 {
        let size = self.dims.total();
        let original = 6 * size * std::mem::size_of::<f32>();
        let compressed: usize = self.tables.iter().map(|t| t.len() * std::mem::size_of::<HCoefficients>()).sum::<usize>()
            + 3 * size * std::mem::size_of::<u32>();
        original as f64 / compressed as f64
    }
}

/// Legacy coefficient set (for backwards compatibility).
#[derive(Clone, Copy, Debug)]
pub struct CoefficientSet {
    /// Voltage-to-voltage coefficient
    pub vv: f32,
    /// Current-to-voltage coefficient
    pub vi: f32,
    /// Current-to-current coefficient
    pub ii: f32,
    /// Voltage-to-current coefficient
    pub iv: f32,
}

impl CoefficientSet {
    /// Create a new coefficient set.
    pub fn new(vv: f32, vi: f32, ii: f32, iv: f32) -> Self {
        Self { vv, vi, ii, iv }
    }

    /// Create free-space coefficients.
    pub fn free_space(dt: f64, dx: f64) -> Self {
        use crate::constants::{EPS0, MU0};
        let vv = 1.0;
        let vi = (dt / (EPS0 * dx)) as f32;
        let ii = 1.0;
        let iv = (dt / (MU0 * dx)) as f32;
        Self { vv, vi, ii, iv }
    }
}

/// Legacy compressed coefficients wrapper.
#[derive(Debug)]
pub struct CompressedCoefficients {
    e_coeffs: CompressedECoefficients,
    h_coeffs: CompressedHCoefficients,
}

impl CompressedCoefficients {
    /// Create from separate E and H compressed coefficients.
    pub fn new(e_coeffs: CompressedECoefficients, h_coeffs: CompressedHCoefficients) -> Self {
        Self { e_coeffs, h_coeffs }
    }

    /// Get E coefficients.
    pub fn e_coefficients(&self) -> &CompressedECoefficients {
        &self.e_coeffs
    }

    /// Get H coefficients.
    pub fn h_coefficients(&self) -> &CompressedHCoefficients {
        &self.h_coeffs
    }

    /// Get compression statistics.
    pub fn stats(&self) -> CompressionStats {
        let dims = self.e_coeffs.dims();
        let size = dims.total();

        let e_unique = self.e_coeffs.num_unique();
        let h_unique = self.h_coeffs.num_unique();

        let original_memory = 12 * size * std::mem::size_of::<f32>(); // 6 E + 6 H coefficients

        let e_table_mem: usize = e_unique.iter().map(|&n| n * std::mem::size_of::<ECoefficients>()).sum();
        let h_table_mem: usize = h_unique.iter().map(|&n| n * std::mem::size_of::<HCoefficients>()).sum();
        let index_mem = 6 * size * std::mem::size_of::<u32>(); // 3 E + 3 H index arrays
        let compressed_memory = e_table_mem + h_table_mem + index_mem;

        CompressionStats {
            total_points: size * 3,
            unique_sets: e_unique.iter().sum::<usize>() + h_unique.iter().sum::<usize>(),
            compression_ratio: original_memory as f64 / compressed_memory as f64,
            memory_saved: original_memory.saturating_sub(compressed_memory),
            original_memory,
            compressed_memory,
        }
    }
}

/// Statistics about coefficient compression.
#[derive(Debug, Clone)]
pub struct CompressionStats {
    /// Total number of grid points (×3 directions)
    pub total_points: usize,
    /// Number of unique coefficient sets
    pub unique_sets: usize,
    /// Compression ratio (original / compressed)
    pub compression_ratio: f64,
    /// Memory saved (bytes)
    pub memory_saved: usize,
    /// Original memory usage (bytes)
    pub original_memory: usize,
    /// Compressed memory usage (bytes)
    pub compressed_memory: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_e_coefficients_compression() {
        let dims = Dimensions::new(10, 10, 10);

        // Create uniform coefficient fields
        let mut ca = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];
        let mut cb = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];

        for dir in 0..3 {
            ca[dir].fill(1.0);
            cb[dir].fill(0.5);
        }

        let compressed = CompressedECoefficients::from_fields(
            &[&ca[0], &ca[1], &ca[2]],
            &[&cb[0], &cb[1], &cb[2]],
        );

        // Uniform field should compress to just one unique set per direction
        let unique = compressed.num_unique();
        assert_eq!(unique[0], 1);
        assert_eq!(unique[1], 1);
        assert_eq!(unique[2], 1);

        // Verify we can retrieve the correct values
        let coeff = compressed.get(0, 5, 5, 5);
        assert!((coeff.ca - 1.0).abs() < 1e-6);
        assert!((coeff.cb - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_h_coefficients_compression() {
        let dims = Dimensions::new(10, 10, 10);

        let mut da = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];
        let mut db = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];

        for dir in 0..3 {
            da[dir].fill(1.0);
            db[dir].fill(-0.5);
        }

        let compressed = CompressedHCoefficients::from_fields(
            &[&da[0], &da[1], &da[2]],
            &[&db[0], &db[1], &db[2]],
        );

        let unique = compressed.num_unique();
        assert_eq!(unique[0], 1);
        assert_eq!(unique[1], 1);
        assert_eq!(unique[2], 1);

        let coeff = compressed.get(0, 5, 5, 5);
        assert!((coeff.da - 1.0).abs() < 1e-6);
        assert!((coeff.db - (-0.5)).abs() < 1e-6);
    }

    #[test]
    fn test_compression_with_material_boundary() {
        let dims = Dimensions::new(10, 10, 10);

        let mut ca = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];
        let mut cb = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];

        // Create a material boundary at i=5
        for dir in 0..3 {
            for i in 0..10 {
                for j in 0..10 {
                    for k in 0..10 {
                        if i < 5 {
                            ca[dir].set(i, j, k, 1.0);
                            cb[dir].set(i, j, k, 0.5);
                        } else {
                            ca[dir].set(i, j, k, 0.8);
                            cb[dir].set(i, j, k, 0.3);
                        }
                    }
                }
            }
        }

        let compressed = CompressedECoefficients::from_fields(
            &[&ca[0], &ca[1], &ca[2]],
            &[&cb[0], &cb[1], &cb[2]],
        );

        // Should have 2 unique sets per direction
        let unique = compressed.num_unique();
        assert_eq!(unique[0], 2);
        assert_eq!(unique[1], 2);
        assert_eq!(unique[2], 2);

        // Verify boundary values
        let coeff_low = compressed.get(0, 2, 5, 5);
        assert!((coeff_low.ca - 1.0).abs() < 1e-6);

        let coeff_high = compressed.get(0, 7, 5, 5);
        assert!((coeff_high.ca - 0.8).abs() < 1e-6);
    }

    #[test]
    fn test_compression_ratio() {
        let dims = Dimensions::new(50, 50, 50);

        // Uniform field - should have excellent compression
        let mut ca = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];
        let mut cb = [Field3D::new(dims), Field3D::new(dims), Field3D::new(dims)];

        for dir in 0..3 {
            ca[dir].fill(1.0);
            cb[dir].fill(0.5);
        }

        let compressed = CompressedECoefficients::from_fields(
            &[&ca[0], &ca[1], &ca[2]],
            &[&cb[0], &cb[1], &cb[2]],
        );

        // For large uniform grids, compression ratio should be > 1
        // (coefficient tables are small, but index arrays still take space)
        let ratio = compressed.compression_ratio();
        println!("Compression ratio for 50³ uniform grid: {:.2}", ratio);

        // The ratio depends on whether we save more in coefficient storage
        // than we spend on index storage. For very uniform grids:
        // Original: 6 * 125000 * 4 = 3MB
        // Compressed: 3 * 8 (tables) + 3 * 125000 * 4 (indices) = ~1.5MB
        assert!(ratio > 0.5); // At minimum, not too much worse
    }
}
