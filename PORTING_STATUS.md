# openEMS C++ to Rust Porting Status

This document tracks the porting progress from C++ to Rust, including implementation notes, file references, and test status.

## Summary

- **Total C++ Extensions**: 11
- **Ported**: 9 (UPML, Mur ABC, Dispersive, Lorentz/Drude/Debye, TF/SF, Lumped RLC, Steady-State, Conducting Sheet, Cylindrical)
- **Remaining**: 2 (Local ABC, Excitation Extension - exists in core)
- **Total Processing Modules**: 7
- **Ported**: 7 (Field probes, Voltage probes, Current probes, VTK output, Integral base, SAR, Mode matching)
- **Remaining**: 0
- **Total Unit Tests**: 95

---

## FDTD Extensions

### 1. UPML (Uniaxial PML) - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_upml.h`
- `FDTD/extensions/operator_ext_upml.cpp`
- `FDTD/extensions/engine_ext_upml.h`
- `FDTD/extensions/engine_ext_upml.cpp`

**Rust Files**:
- `src/extensions/pml.rs`

**Implementation Notes**:
- Flux-based ADE formulation from Taflove Chapter 7.8
- Polynomial grading with configurable order
- Per-boundary configuration (enable/disable each of 6 faces)
- Pre/post voltage and current update hooks
- Coefficient arrays: vv, vvfo, vvfn, ii, iifo, iifn

**Tests**:
- `test_pml_config` - Configuration parameters
- `test_pml_boundaries` - Boundary enable/disable
- `test_upml_creation` - Region initialization
- `test_upml_disabled_boundary` - Partial PML
- `test_kappa_grading` - Conductivity profile
- `test_upml_field_update` - Update cycle

---

### 2. Mur ABC (1st Order Absorbing BC) - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_mur_abc.h`
- `FDTD/extensions/operator_ext_mur_abc.cpp`
- `FDTD/extensions/engine_ext_mur_abc.h`
- `FDTD/extensions/engine_ext_mur_abc.cpp`

**Rust Files**:
- `src/extensions/mur_abc.rs`

**Implementation Notes**:
- First-order Mur ABC based on one-way wave equation
- Per-face coefficient arrays
- Pre/post voltage update hooks
- Formula: Mur_coeff = (c*dt - delta) / (c*dt + delta)
- Configurable phase velocity

**Tests**:
- `test_mur_config` - Configuration
- `test_mur_creation` - Face initialization
- `test_mur_coefficient` - Coefficient calculation
- `test_mur_partial_boundaries` - Selective faces
- `test_mur_update_cycle` - Update sequence
- `test_mur_with_custom_phase_velocity` - Phase velocity override

---

### 3. Local Absorbing BC - PENDING
**Status**: Not started (lower priority - Mur ABC covers most use cases)

**C++ Files**:
- `FDTD/extensions/operator_ext_absorbing_bc.h`
- `FDTD/extensions/operator_ext_absorbing_bc.cpp`

---

### 4. Dispersive Materials (Lorentz/Drude/Debye) - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_dispersive.h`
- `FDTD/extensions/operator_ext_lorentzmaterial.h`
- `FDTD/extensions/engine_ext_lorentzmaterial.cpp`

**Rust Files**:
- `src/extensions/dispersive.rs`

**Implementation Notes**:
- ADE (Auxiliary Differential Equation) formulation
- Three material models:
  - **Lorentz**: Dielectric resonance (oscillator)
  - **Drude**: Metals and plasmas (free electron)
  - **Debye**: Polar molecules (relaxation)
- Position-based material assignment
- Auxiliary current storage per position

**Tests**:
- `test_lorentz_params` - Parameter creation
- `test_drude_params` - Metal parameters
- `test_lorentz_material` - Material update
- `test_drude_material` - Drude update
- `test_debye_material` - Debye update
- `test_material_reset` - State reset

---

### 5. Conducting Sheet Model - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_conductingsheet.h`
- `FDTD/extensions/operator_ext_conductingsheet.cpp`
- `FDTD/extensions/cond_sheet_parameter.h`

**Rust Files**:
- `src/extensions/conducting_sheet.rs`

**Implementation Notes**:
- Based on Lauer & Wolff (1999) "A conducting sheet model for efficient wide band FDTD analysis"
- Frequency-dependent sheet impedance using Lorentz pole decomposition
- ADE coefficients for time-domain update
- Support for copper and aluminum sheets
- Sheet thickness and conductivity parameters
- Manager for multiple sheets

**Tests**:
- `test_lorentz_params` - Lorentz parameter selection
- `test_sheet_config` - Configuration for copper/aluminum
- `test_sheet_creation` - Sheet initialization
- `test_sheet_update` - E-field modification
- `test_sheet_reset` - State reset
- `test_sheet_manager` - Multiple sheet management

---

### 6. Excitation Extension - EXISTS IN CORE
**Status**: Excitation handled in `src/fdtd/excitation.rs`

---

### 7. TF/SF (Total-Field/Scattered-Field) - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_tfsf.h`
- `FDTD/extensions/operator_ext_tfsf.cpp`
- `FDTD/extensions/engine_ext_tfsf.h`
- `FDTD/extensions/engine_ext_tfsf.cpp`

**Rust Files**:
- `src/extensions/tfsf.rs`

**Implementation Notes**:
- 1D auxiliary FDTD grid for plane wave propagation
- Incident field injection at TF/SF boundary
- Configurable propagation direction (6 options)
- Configurable polarization
- ABC at 1D grid boundaries

**Tests**:
- `test_tfsf_config` - Configuration
- `test_propagation_direction` - Direction handling
- `test_tfsf_creation` - Boundary initialization
- `test_incident_field` - Field calculation
- `test_tfsf_update` - Update cycle
- `test_tfsf_reset` - State reset

---

### 8. Lumped RLC Elements - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_lumpedRLC.h`
- `FDTD/extensions/operator_ext_lumpedRLC.cpp`
- `FDTD/extensions/engine_ext_lumpedRLC.h`
- `FDTD/extensions/engine_ext_lumpedRLC.cpp`

**Rust Files**:
- `src/extensions/lumped_rlc.rs`

**Implementation Notes**:
- Pure R, L, C elements
- Parallel and Series RLC combinations
- Update coefficients for E-field modification
- Integral storage for inductance

**Tests**:
- `test_rlc_element_creation` - Element types
- `test_lumped_rlc_creation` - Extension creation
- `test_rlc_update` - Update cycle
- `test_parallel_rlc` - Parallel config
- `test_series_rlc` - Series config
- `test_rlc_reset` - State reset

---

### 9. Cylindrical Coordinate Extension - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/operator_cylinder.h`
- `FDTD/operator_cylinder.cpp`
- `FDTD/engine_cylinder.h`
- `FDTD/engine_cylinder.cpp`
- `FDTD/extensions/operator_ext_cylinder.h`

**Rust Files**:
- `src/fdtd/cylindrical.rs`

**Implementation Notes**:
- Full cylindrical FDTD (rho, alpha, z coordinates)
- Yee grid adapted for cylindrical geometry
- Closed alpha mesh support (full 360 degrees)
- R=0 singularity handling with special coefficients
- Radius-dependent edge lengths and cell areas
- CFL timestep calculation for cylindrical grids

**Tests**:
- `test_cylindrical_grid_creation` - Grid initialization
- `test_grid_with_r0` - R=0 included grids
- `test_cell_volume` - Volume calculations
- `test_alpha_wrapping` - Closed mesh index wrapping
- `test_cylindrical_operator` - Coefficient calculation
- `test_cylindrical_engine_creation` - Engine initialization
- `test_engine_step` - Time stepping
- `test_engine_reset` - State reset
- `test_edge_length_with_radius` - Arc length scaling
- `test_node_area` - Face area calculations

---

### 10. Steady-State Detection - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `FDTD/extensions/operator_ext_steadystate.h`
- `FDTD/extensions/engine_ext_steadystate.cpp`

**Rust Files**:
- `src/extensions/steady_state.rs`

**Implementation Notes**:
- Field probe monitoring at strategic locations
- Period detection via zero-crossing analysis
- Convergence threshold comparison across periods
- Configurable check intervals and minimum timesteps

**Tests**:
- `test_config` - Configuration
- `test_detector_creation` - Detector initialization
- `test_add_probes` - Probe management
- `test_default_probes` - Default placement
- `test_record` - Recording
- `test_not_converged_early` - Early termination prevention
- `test_reset` - State reset

---

## Processing Modules

### 1. Voltage Probe - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `Common/processvoltage.h`
- `Common/processvoltage.cpp`

**Rust Files**:
- `src/processing/probes.rs`

**Implementation Notes**:
- Line integral of E-field
- FFT for frequency response

**Tests**:
- `test_voltage_probe` - Basic functionality

---

### 2. Current Probe - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `Common/processcurrent.h`
- `Common/processcurrent.cpp`

**Rust Files**:
- `src/processing/probes.rs`

**Implementation Notes**:
- Surface integral of H-field
- FFT for frequency response

**Tests**:
- `test_current_probe` - Basic functionality

---

### 3. Field Probes - COMPLETED
**Status**: Implemented with tests

**Rust Files**:
- `src/processing/probes.rs`

**Tests**:
- `test_field_probe` - Point measurement
- `test_processing_collection` - Collection management

---

### 4. Field Dumps (TD/FD) - PARTIAL
**Status**: VTK output exists, HDF5 pending

**Rust Files**:
- `src/io/vtk.rs`

---

### 5. SAR Calculation - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `Common/processfields_sar.h`
- `Common/processfields_sar.cpp`
- `tools/sar_calculation.h`
- `tools/sar_calculation.cpp`

**Rust Files**:
- `src/processing/sar.rs`

**Implementation Notes**:
- Frequency domain field storage (E and J fields)
- Local SAR: P = 0.5 * sigma * |E|^2 / density
- Multiple averaging methods: IEEE C95.3, IEEE 62704, Simple, Local
- Cubical mass fitting for averaged SAR (1g, 10g)
- Total absorbed power calculation
- Statistics tracking (peak SAR, air voxels, etc.)

**Tests**:
- `test_sar_config` - Configuration options
- `test_averaging_method_parse` - Method string parsing
- `test_fd_fields_creation` - Frequency domain storage
- `test_sar_calculation_creation` - Calculator initialization
- `test_local_sar_calculation` - Point SAR values
- `test_averaged_sar` - Mass-averaged SAR
- `test_total_power` - Power integration
- `test_sar_stats` - Statistics tracking

---

### 6. Mode Matching - COMPLETED
**Status**: Implemented with tests

**C++ Files**:
- `Common/processmodematch.h`
- `Common/processmodematch.cpp`

**Rust Files**:
- `src/processing/mode_match.rs`

**Implementation Notes**:
- Mode function trait for arbitrary mode definitions
- Built-in rectangular waveguide TE/TM modes
- Built-in circular waveguide TE modes
- Custom mode function support via closures
- Mode integral (voltage/current) calculation
- Mode purity calculation
- Time series recording

**Tests**:
- `test_rectangular_te_mode` - TE mode functions
- `test_mode_match_config` - Configuration
- `test_mode_match_creation` - Matcher initialization
- `test_mode_match_with_function` - Full mode matching
- `test_custom_mode_function` - User-defined modes
- `test_mode_match_results` - Time series recording
- `test_circular_te_mode` - Circular waveguide modes

---

## Engine Implementations

### 1. Basic Cartesian - COMPLETED
**Rust**: `src/fdtd/engine.rs` (Basic variant)

### 2. SIMD-Optimized - COMPLETED
**Rust**: `src/fdtd/engine.rs` (Simd variant)

### 3. Multi-threaded - COMPLETED
**Rust**: `src/fdtd/engine.rs` (Parallel variant using rayon)

### 4. Cylindrical Coordinates - COMPLETED
**Rust**: `src/fdtd/cylindrical.rs`

### 5. MPI Distributed - NOT PLANNED
**Note**: Out of scope for Rust port (use alternative parallelism)

---

## Core Modules

### 1. Operator (FDTD Coefficients) - COMPLETED
**Rust**: `src/fdtd/operator.rs`

### 2. Excitation (Sources) - COMPLETED
**Rust**: `src/fdtd/excitation.rs`

### 3. Grid/Geometry - COMPLETED
**Rust**: `src/geometry/mod.rs`

### 4. Arrays/Fields - COMPLETED
**Rust**: `src/arrays/mod.rs`, `src/arrays/field.rs`

### 5. Signal Processing - COMPLETED
**Rust**: `src/tools/signal.rs`

### 6. NF2FF Transform - COMPLETED
**Rust**: `src/nf2ff/mod.rs`

---

## I/O Modules

### 1. VTK Output - COMPLETED
**Rust**: `src/io/vtk.rs`

### 2. XML Input - COMPLETED
**Rust**: `src/io/xml.rs`

### 3. HDF5 I/O - PENDING
**C++**: `tools/hdf5_file_writer.h`, `tools/hdf5_file_reader.h`

---

## Changelog

### 2026-01-25 (Session 2)
- Added SAR calculation processor with IEEE averaging methods
- Added Mode Matching processor with rectangular/circular waveguide modes
- Added Conducting Sheet extension with frequency-dependent impedance
- Added Cylindrical Coordinate FDTD implementation
- Added Debug trait and zero() method to VectorField3D
- Total tests increased from 65 to 95

### 2026-01-25 (Session 1)
- Created tracking document
- Completed UPML implementation with full test coverage
- Completed Mur ABC implementation with tests
- Completed dispersive materials (Lorentz/Drude/Debye) with tests
- Completed TF/SF boundary with tests
- Completed Lumped RLC elements with tests
- Completed Steady-State Detection with tests
- Completed Voltage and Current probes with tests
- Added 65 passing unit tests total
