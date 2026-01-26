//! HDF5 file I/O support.
//!
//! Provides reading and writing of simulation data in HDF5 format.
//! This module requires the `hdf5` feature to be enabled.
//!
//! When HDF5 is not available, a fallback using memory-mapped binary files is provided.

use crate::arrays::{Dimensions, VectorField3D};
use num_complex::Complex64;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// Mesh type enumeration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeshType {
    /// Cartesian mesh
    Cartesian = 0,
    /// Cylindrical mesh
    Cylindrical = 1,
}

impl From<i32> for MeshType {
    fn from(v: i32) -> Self {
        match v {
            1 => MeshType::Cylindrical,
            _ => MeshType::Cartesian,
        }
    }
}

/// HDF5 file reader for openEMS data.
#[derive(Debug)]
#[allow(dead_code)]
pub struct Hdf5Reader {
    /// File path
    path: String,
    /// Mesh lines for each direction
    mesh_lines: [Vec<f64>; 3],
    /// Mesh type
    mesh_type: MeshType,
    /// Timesteps (for TD data)
    timesteps: Vec<f64>,
    /// Frequencies (for FD data)
    frequencies: Vec<f64>,
    /// Attributes
    attributes: HashMap<String, Vec<f64>>,
}

impl Hdf5Reader {
    /// Open an HDF5 file for reading.
    pub fn open<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let path_str = path.as_ref().to_string_lossy().to_string();

        // Try to read as our binary format first
        if let Ok(reader) = Self::read_binary_format(&path_str) {
            return Ok(reader);
        }

        #[cfg(feature = "hdf5")]
        {
            return Self::read_hdf5_format(&path_str);
        }

        #[cfg(not(feature = "hdf5"))]
        {
            Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "HDF5 support not enabled. Enable the 'hdf5' feature or use binary format.",
            ))
        }
    }

    /// Read our custom binary format.
    fn read_binary_format(path: &str) -> io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read magic number
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if &magic != b"OEMS" {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "Invalid file format"));
        }

        // Read version
        let mut version = [0u8; 4];
        reader.read_exact(&mut version)?;

        // Read mesh type
        let mut mesh_type_buf = [0u8; 4];
        reader.read_exact(&mut mesh_type_buf)?;
        let mesh_type = MeshType::from(i32::from_le_bytes(mesh_type_buf));

        // Read mesh lines
        let mut mesh_lines = [Vec::new(), Vec::new(), Vec::new()];
        for lines in &mut mesh_lines {
            let mut len_buf = [0u8; 8];
            reader.read_exact(&mut len_buf)?;
            let len = u64::from_le_bytes(len_buf) as usize;

            *lines = vec![0.0; len];
            for val in lines.iter_mut() {
                let mut val_buf = [0u8; 8];
                reader.read_exact(&mut val_buf)?;
                *val = f64::from_le_bytes(val_buf);
            }
        }

        // Read timesteps
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let ts_len = u64::from_le_bytes(len_buf) as usize;
        let mut timesteps = vec![0.0; ts_len];
        for ts in &mut timesteps {
            let mut val_buf = [0u8; 8];
            reader.read_exact(&mut val_buf)?;
            *ts = f64::from_le_bytes(val_buf);
        }

        // Read frequencies
        let mut len_buf = [0u8; 8];
        reader.read_exact(&mut len_buf)?;
        let freq_len = u64::from_le_bytes(len_buf) as usize;
        let mut frequencies = vec![0.0; freq_len];
        for freq in &mut frequencies {
            let mut val_buf = [0u8; 8];
            reader.read_exact(&mut val_buf)?;
            *freq = f64::from_le_bytes(val_buf);
        }

        Ok(Self {
            path: path.to_string(),
            mesh_lines,
            mesh_type,
            timesteps,
            frequencies,
            attributes: HashMap::new(),
        })
    }

    #[cfg(feature = "hdf5")]
    fn read_hdf5_format(path: &str) -> io::Result<Self> {
        use hdf5_metno::File as H5File;

        let file = H5File::open(path).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
        })?;

        // Read mesh
        let mut mesh_lines = [Vec::new(), Vec::new(), Vec::new()];
        if let Ok(mesh_group) = file.group("Mesh") {
            for (i, name) in ["x", "y", "z"].iter().enumerate() {
                if let Ok(dataset) = mesh_group.dataset(name) {
                    if let Ok(data) = dataset.read_1d::<f64>() {
                        mesh_lines[i] = data.to_vec();
                    }
                }
            }
        }

        // Determine mesh type from attribute
        let mesh_type = if let Ok(attr) = file.attr("MeshType") {
            if let Ok(val) = attr.read_scalar::<i32>() {
                MeshType::from(val)
            } else {
                MeshType::Cartesian
            }
        } else {
            MeshType::Cartesian
        };

        // Read timesteps
        let timesteps = if let Ok(td_group) = file.group("FieldData/TD") {
            if let Ok(dataset) = td_group.dataset("time") {
                dataset.read_1d::<f64>().map(|d| d.to_vec()).unwrap_or_default()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        // Read frequencies
        let frequencies = if let Ok(fd_group) = file.group("FieldData/FD") {
            if let Ok(dataset) = fd_group.dataset("frequency") {
                dataset.read_1d::<f64>().map(|d| d.to_vec()).unwrap_or_default()
            } else {
                Vec::new()
            }
        } else {
            Vec::new()
        };

        Ok(Self {
            path: path.to_string(),
            mesh_lines,
            mesh_type,
            timesteps,
            frequencies,
            attributes: HashMap::new(),
        })
    }

    /// Get mesh lines.
    pub fn mesh_lines(&self) -> &[Vec<f64>; 3] {
        &self.mesh_lines
    }

    /// Get mesh type.
    pub fn mesh_type(&self) -> MeshType {
        self.mesh_type
    }

    /// Get number of timesteps.
    pub fn num_timesteps(&self) -> usize {
        self.timesteps.len()
    }

    /// Get timesteps.
    pub fn timesteps(&self) -> &[f64] {
        &self.timesteps
    }

    /// Get number of frequencies.
    pub fn num_frequencies(&self) -> usize {
        self.frequencies.len()
    }

    /// Get frequencies.
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }

    /// Check if file is valid.
    pub fn is_valid(&self) -> bool {
        !self.mesh_lines[0].is_empty()
            && !self.mesh_lines[1].is_empty()
            && !self.mesh_lines[2].is_empty()
    }

    /// Get dimensions from mesh.
    pub fn dimensions(&self) -> Dimensions {
        Dimensions {
            nx: self.mesh_lines[0].len().saturating_sub(1).max(1),
            ny: self.mesh_lines[1].len().saturating_sub(1).max(1),
            nz: self.mesh_lines[2].len().saturating_sub(1).max(1),
        }
    }
}

/// HDF5 file writer for openEMS data.
#[derive(Debug)]
pub struct Hdf5Writer {
    /// File path
    path: String,
    /// Current group path
    current_group: String,
    /// Mesh lines
    mesh_lines: Option<[Vec<f64>; 3]>,
    /// Mesh type
    mesh_type: MeshType,
    /// Pending scalar fields
    scalar_fields: Vec<(String, Vec<f32>)>,
    /// Pending vector fields
    vector_fields: Vec<(String, [Vec<f32>; 3])>,
    /// Pending complex fields
    complex_fields: Vec<(String, Vec<Complex64>)>,
    /// Attributes
    attributes: HashMap<String, Vec<f64>>,
}

impl Hdf5Writer {
    /// Create a new HDF5 writer.
    pub fn new<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Self {
            path: path.as_ref().to_string_lossy().to_string(),
            current_group: "/".to_string(),
            mesh_lines: None,
            mesh_type: MeshType::Cartesian,
            scalar_fields: Vec::new(),
            vector_fields: Vec::new(),
            complex_fields: Vec::new(),
            attributes: HashMap::new(),
        })
    }

    /// Set the current group for writing.
    pub fn set_current_group(&mut self, group: &str) {
        self.current_group = group.to_string();
    }

    /// Write rectangular mesh.
    pub fn write_rect_mesh(&mut self, lines: [Vec<f64>; 3], mesh_type: MeshType) {
        self.mesh_lines = Some(lines);
        self.mesh_type = mesh_type;
    }

    /// Write a scalar field.
    pub fn write_scalar_field(&mut self, name: &str, data: &[f32]) {
        self.scalar_fields.push((name.to_string(), data.to_vec()));
    }

    /// Write a scalar field (f64).
    pub fn write_scalar_field_f64(&mut self, name: &str, data: &[f64]) {
        let f32_data: Vec<f32> = data.iter().map(|&x| x as f32).collect();
        self.scalar_fields.push((name.to_string(), f32_data));
    }

    /// Write a vector field.
    pub fn write_vector_field(&mut self, name: &str, field: &VectorField3D) {
        let x_data = field.x.as_slice().to_vec();
        let y_data = field.y.as_slice().to_vec();
        let z_data = field.z.as_slice().to_vec();
        self.vector_fields.push((name.to_string(), [x_data, y_data, z_data]));
    }

    /// Write a complex scalar field.
    pub fn write_complex_field(&mut self, name: &str, data: &[Complex64]) {
        self.complex_fields.push((name.to_string(), data.to_vec()));
    }

    /// Write an attribute.
    pub fn write_attribute(&mut self, name: &str, values: &[f64]) {
        self.attributes.insert(name.to_string(), values.to_vec());
    }

    /// Flush and close the file.
    pub fn close(self) -> io::Result<()> {
        #[cfg(feature = "hdf5")]
        {
            return self.write_hdf5();
        }

        #[cfg(not(feature = "hdf5"))]
        {
            self.write_binary()
        }
    }

    /// Write in our custom binary format.
    fn write_binary(self) -> io::Result<()> {
        let file = File::create(&self.path)?;
        let mut writer = BufWriter::new(file);

        // Write magic number
        writer.write_all(b"OEMS")?;

        // Write version (1.0)
        writer.write_all(&[1u8, 0, 0, 0])?;

        // Write mesh type
        writer.write_all(&(self.mesh_type as i32).to_le_bytes())?;

        // Write mesh lines
        if let Some(ref mesh_lines) = self.mesh_lines {
            for lines in mesh_lines.iter() {
                writer.write_all(&(lines.len() as u64).to_le_bytes())?;
                for &val in lines.iter() {
                    writer.write_all(&val.to_le_bytes())?;
                }
            }
        } else {
            // Write empty mesh
            for _ in 0..3 {
                writer.write_all(&0u64.to_le_bytes())?;
            }
        }

        // Write timesteps placeholder
        writer.write_all(&0u64.to_le_bytes())?;

        // Write frequencies placeholder
        writer.write_all(&0u64.to_le_bytes())?;

        // Write scalar fields
        writer.write_all(&(self.scalar_fields.len() as u64).to_le_bytes())?;
        for (name, data) in &self.scalar_fields {
            // Write name length and name
            writer.write_all(&(name.len() as u64).to_le_bytes())?;
            writer.write_all(name.as_bytes())?;
            // Write data
            writer.write_all(&(data.len() as u64).to_le_bytes())?;
            for &val in data.iter() {
                writer.write_all(&val.to_le_bytes())?;
            }
        }

        // Write vector fields
        writer.write_all(&(self.vector_fields.len() as u64).to_le_bytes())?;
        for (name, components) in &self.vector_fields {
            writer.write_all(&(name.len() as u64).to_le_bytes())?;
            writer.write_all(name.as_bytes())?;
            for component in components.iter() {
                writer.write_all(&(component.len() as u64).to_le_bytes())?;
                for &val in component.iter() {
                    writer.write_all(&val.to_le_bytes())?;
                }
            }
        }

        // Write complex fields
        writer.write_all(&(self.complex_fields.len() as u64).to_le_bytes())?;
        for (name, data) in &self.complex_fields {
            writer.write_all(&(name.len() as u64).to_le_bytes())?;
            writer.write_all(name.as_bytes())?;
            writer.write_all(&(data.len() as u64).to_le_bytes())?;
            for &val in data.iter() {
                writer.write_all(&val.re.to_le_bytes())?;
                writer.write_all(&val.im.to_le_bytes())?;
            }
        }

        writer.flush()?;
        Ok(())
    }

    #[cfg(feature = "hdf5")]
    fn write_hdf5(self) -> io::Result<()> {
        use hdf5_metno::File as H5File;

        let file = H5File::create(&self.path).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
        })?;

        // Write mesh
        if let Some(ref mesh_lines) = self.mesh_lines {
            let mesh_group = file.create_group("Mesh").map_err(|e| {
                io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
            })?;

            for (i, name) in ["x", "y", "z"].iter().enumerate() {
                let dataset = mesh_group.new_dataset::<f64>()
                    .shape([mesh_lines[i].len()])
                    .create(name)
                    .map_err(|e| {
                        io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
                    })?;
                dataset.write(&mesh_lines[i]).map_err(|e| {
                    io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
                })?;
            }
        }

        // Write mesh type attribute
        let attr = file.new_attr::<i32>().create("MeshType").map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
        })?;
        attr.write_scalar(&(self.mesh_type as i32)).map_err(|e| {
            io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
        })?;

        // Write scalar fields
        for (name, data) in &self.scalar_fields {
            let dataset = file.new_dataset::<f32>()
                .shape([data.len()])
                .create(&name)
                .map_err(|e| {
                    io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
                })?;
            dataset.write(data).map_err(|e| {
                io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
            })?;
        }

        // Write vector fields
        for (name, components) in &self.vector_fields {
            let group = file.create_group(&name).map_err(|e| {
                io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
            })?;

            for (i, comp_name) in ["x", "y", "z"].iter().enumerate() {
                let dataset = group.new_dataset::<f32>()
                    .shape([components[i].len()])
                    .create(comp_name)
                    .map_err(|e| {
                        io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
                    })?;
                dataset.write(&components[i]).map_err(|e| {
                    io::Error::new(io::ErrorKind::Other, format!("HDF5 error: {}", e))
                })?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn test_binary_roundtrip() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.oems");

        // Write
        let mut writer = Hdf5Writer::new(&path).unwrap();
        writer.write_rect_mesh(
            [
                vec![0.0, 0.001, 0.002],
                vec![0.0, 0.001, 0.002],
                vec![0.0, 0.001, 0.002],
            ],
            MeshType::Cartesian,
        );
        writer.write_scalar_field("test_field", &[1.0, 2.0, 3.0, 4.0]);
        writer.close().unwrap();

        // Read
        let reader = Hdf5Reader::open(&path).unwrap();
        assert!(reader.is_valid());
        assert_eq!(reader.mesh_type(), MeshType::Cartesian);
        assert_eq!(reader.mesh_lines()[0].len(), 3);
    }

    #[test]
    fn test_mesh_type_conversion() {
        assert_eq!(MeshType::from(0), MeshType::Cartesian);
        assert_eq!(MeshType::from(1), MeshType::Cylindrical);
        assert_eq!(MeshType::from(99), MeshType::Cartesian);
    }

    #[test]
    fn test_writer_creation() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test.h5");

        let writer = Hdf5Writer::new(&path).unwrap();
        assert_eq!(writer.current_group, "/");
    }
}
