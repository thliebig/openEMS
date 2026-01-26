//! VTK file format support for visualization.

use crate::arrays::Field3D;
use crate::geometry::Grid;
use crate::Result;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Write a 3D field to VTK format (legacy ASCII).
pub fn write_field_vtk(
    path: &Path,
    field: &Field3D,
    grid: &Grid,
    field_name: &str,
) -> Result<()> {
    let dims = field.dims();
    let mut file = File::create(path)?;

    // VTK header
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "openEMS field data")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET RECTILINEAR_GRID")?;
    writeln!(
        file,
        "DIMENSIONS {} {} {}",
        dims.nx + 1,
        dims.ny + 1,
        dims.nz + 1
    )?;

    // X coordinates
    writeln!(file, "X_COORDINATES {} float", dims.nx + 1)?;
    for x in grid.x_lines() {
        write!(file, "{} ", x)?;
    }
    writeln!(file)?;

    // Y coordinates
    writeln!(file, "Y_COORDINATES {} float", dims.ny + 1)?;
    for y in grid.y_lines() {
        write!(file, "{} ", y)?;
    }
    writeln!(file)?;

    // Z coordinates
    writeln!(file, "Z_COORDINATES {} float", dims.nz + 1)?;
    for z in grid.z_lines() {
        write!(file, "{} ", z)?;
    }
    writeln!(file)?;

    // Cell data
    writeln!(file, "CELL_DATA {}", dims.total())?;
    writeln!(file, "SCALARS {} float 1", field_name)?;
    writeln!(file, "LOOKUP_TABLE default")?;

    for i in 0..dims.nx {
        for j in 0..dims.ny {
            for k in 0..dims.nz {
                writeln!(file, "{}", field.get(i, j, k))?;
            }
        }
    }

    Ok(())
}

/// Write a vector field to VTK format.
pub fn write_vector_field_vtk(
    path: &Path,
    field: &crate::arrays::VectorField3D,
    grid: &Grid,
    field_name: &str,
) -> Result<()> {
    let dims = field.dims();
    let mut file = File::create(path)?;

    // VTK header
    writeln!(file, "# vtk DataFile Version 3.0")?;
    writeln!(file, "openEMS vector field data")?;
    writeln!(file, "ASCII")?;
    writeln!(file, "DATASET RECTILINEAR_GRID")?;
    writeln!(
        file,
        "DIMENSIONS {} {} {}",
        dims.nx + 1,
        dims.ny + 1,
        dims.nz + 1
    )?;

    // Coordinates
    writeln!(file, "X_COORDINATES {} float", dims.nx + 1)?;
    for x in grid.x_lines() {
        write!(file, "{} ", x)?;
    }
    writeln!(file)?;

    writeln!(file, "Y_COORDINATES {} float", dims.ny + 1)?;
    for y in grid.y_lines() {
        write!(file, "{} ", y)?;
    }
    writeln!(file)?;

    writeln!(file, "Z_COORDINATES {} float", dims.nz + 1)?;
    for z in grid.z_lines() {
        write!(file, "{} ", z)?;
    }
    writeln!(file)?;

    // Vector data
    writeln!(file, "CELL_DATA {}", dims.total())?;
    writeln!(file, "VECTORS {} float", field_name)?;

    for i in 0..dims.nx {
        for j in 0..dims.ny {
            for k in 0..dims.nz {
                writeln!(
                    file,
                    "{} {} {}",
                    field.x.get(i, j, k),
                    field.y.get(i, j, k),
                    field.z.get(i, j, k)
                )?;
            }
        }
    }

    Ok(())
}
