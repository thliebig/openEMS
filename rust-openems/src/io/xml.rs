//! XML configuration file support.

use crate::fdtd::BoundaryCondition;
use crate::{Error, Result};
use quick_xml::events::{BytesStart, Event};
use quick_xml::Reader;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Parse an openEMS XML configuration file.
pub fn parse_config(path: &Path) -> Result<SimulationConfig> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut xml_reader = Reader::from_reader(reader);
    xml_reader.config_mut().trim_text(true);

    let mut config = SimulationConfig::default();
    let mut buf = Vec::new();

    loop {
        match xml_reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                match e.name().as_ref() {
                    b"openEMS" => {
                        // Root element
                    }
                    b"FDTD" => {
                        parse_fdtd_settings(&mut xml_reader, e, &mut config)?;
                    }
                    b"ContinuousStructure" => {
                        parse_geometry(&mut xml_reader, &mut config)?;
                    }
                    _ => {}
                }
            }
            Ok(Event::Eof) => break,
            Err(e) => return Err(Error::Xml(format!("XML parse error: {}", e))),
            _ => {}
        }
        buf.clear();
    }

    Ok(config)
}

/// Simulation configuration parsed from XML.
#[derive(Debug, Default)]
pub struct SimulationConfig {
    /// Maximum number of timesteps
    pub num_timesteps: u64,
    /// End criteria (energy decay in dB)
    pub end_criteria_db: Option<f64>,
    /// Excitation type
    pub excitation_type: ExcitationType,
    /// Center frequency for excitation
    pub f0: f64,
    /// Cutoff frequency
    pub fc: f64,
    /// Boundary conditions
    pub boundaries: [i32; 6],
    /// Grid lines
    pub grid: Option<GridConfig>,
}

/// Excitation configuration.
#[derive(Debug, Default)]
pub enum ExcitationType {
    #[default]
    Gaussian,
    Sinusoidal,
    Dirac,
    Step,
    Custom,
}

/// Grid configuration.
#[derive(Debug, Default)]
pub struct GridConfig {
    pub x_lines: Vec<f64>,
    pub y_lines: Vec<f64>,
    pub z_lines: Vec<f64>,
    pub unit: f64,
}

fn parse_fdtd_settings(
    reader: &mut Reader<BufReader<File>>,
    start: &BytesStart,
    config: &mut SimulationConfig,
) -> Result<()> {
    // Parse FDTD attributes
    for attr in start.attributes() {
        if let Ok(attr) = attr {
            match attr.key.as_ref() {
                b"NumberOfTimesteps" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.num_timesteps = val.parse().unwrap_or(10000);
                    }
                }
                b"endCriteria" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.end_criteria_db = val.parse().ok();
                    }
                }
                b"f0" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.f0 = val.parse().unwrap_or(0.0);
                    }
                }
                b"fc" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.fc = val.parse().unwrap_or(0.0);
                    }
                }
                _ => {}
            }
        }
    }

    // Parse child elements
    let mut buf = Vec::new();
    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(ref e)) => {
                if e.name().as_ref() == b"BoundaryCond" {
                    parse_boundary_cond(e, config)?;
                }
            }
            Ok(Event::End(ref e)) if e.name().as_ref() == b"FDTD" => break,
            Ok(Event::Eof) => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(())
}

fn parse_boundary_cond(start: &BytesStart, config: &mut SimulationConfig) -> Result<()> {
    for attr in start.attributes() {
        if let Ok(attr) = attr {
            match attr.key.as_ref() {
                b"xmin" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.boundaries[0] = val.parse().unwrap_or(0);
                    }
                }
                b"xmax" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.boundaries[1] = val.parse().unwrap_or(0);
                    }
                }
                b"ymin" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.boundaries[2] = val.parse().unwrap_or(0);
                    }
                }
                b"ymax" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.boundaries[3] = val.parse().unwrap_or(0);
                    }
                }
                b"zmin" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.boundaries[4] = val.parse().unwrap_or(0);
                    }
                }
                b"zmax" => {
                    if let Ok(val) = std::str::from_utf8(&attr.value) {
                        config.boundaries[5] = val.parse().unwrap_or(0);
                    }
                }
                _ => {}
            }
        }
    }
    Ok(())
}

fn parse_geometry(
    reader: &mut Reader<BufReader<File>>,
    _config: &mut SimulationConfig,
) -> Result<()> {
    // TODO: Parse geometry elements (primitives, materials, etc.)
    let mut buf = Vec::new();
    let mut depth = 1;

    loop {
        match reader.read_event_into(&mut buf) {
            Ok(Event::Start(_)) => depth += 1,
            Ok(Event::End(_)) => {
                depth -= 1;
                if depth == 0 {
                    break;
                }
            }
            Ok(Event::Eof) => break,
            _ => {}
        }
        buf.clear();
    }

    Ok(())
}

/// Convert boundary condition code to enum.
pub fn boundary_from_code(code: i32) -> BoundaryCondition {
    match code {
        0 => BoundaryCondition::Pec,
        1 => BoundaryCondition::Pmc,
        2 => BoundaryCondition::MurAbc,
        3 => BoundaryCondition::Pml { layers: 8 },
        4 => BoundaryCondition::Periodic,
        _ => BoundaryCondition::Pec,
    }
}
