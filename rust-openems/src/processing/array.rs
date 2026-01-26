//! Processing array - container for multiple processing objects.
//!
//! Manages a collection of processing objects and handles their lifecycle.

use crate::arrays::VectorField3D;
use std::any::Any;

/// Trait for all processing objects.
pub trait Processor: Send + Sync {
    /// Get the processor name.
    fn name(&self) -> &str;

    /// Initialize the processor.
    fn init(&mut self) {}

    /// Pre-process step (called before simulation starts).
    fn pre_process(&mut self) {}

    /// Process a timestep.
    fn process(
        &mut self,
        e_field: &VectorField3D,
        h_field: &VectorField3D,
        timestep: u64,
        time: f64,
    );

    /// Post-process step (called after simulation ends).
    fn post_process(&mut self) {}

    /// Flush data to disk.
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }

    /// Check if processing is needed at this timestep.
    fn needs_processing(&self, _timestep: u64) -> bool {
        true
    }

    /// Get processing interval (0 = every timestep).
    fn interval(&self) -> u64 {
        0
    }

    /// Reset the processor state.
    fn reset(&mut self) {}

    /// As any for downcasting.
    fn as_any(&self) -> &dyn Any;

    /// As any mut for downcasting.
    fn as_any_mut(&mut self) -> &mut dyn Any;
}

/// Container for multiple processing objects.
pub struct ProcessingArray {
    /// Collection of processors
    processors: Vec<Box<dyn Processor>>,
    /// Whether the array has been initialized
    initialized: bool,
}

impl ProcessingArray {
    /// Create a new empty processing array.
    pub fn new() -> Self {
        Self {
            processors: Vec::new(),
            initialized: false,
        }
    }

    /// Add a processor to the array.
    pub fn add<P: Processor + 'static>(&mut self, processor: P) {
        self.processors.push(Box::new(processor));
    }

    /// Add a boxed processor to the array.
    pub fn add_boxed(&mut self, processor: Box<dyn Processor>) {
        self.processors.push(processor);
    }

    /// Get number of processors.
    pub fn len(&self) -> usize {
        self.processors.len()
    }

    /// Check if array is empty.
    pub fn is_empty(&self) -> bool {
        self.processors.is_empty()
    }

    /// Initialize all processors.
    pub fn init_all(&mut self) {
        for proc in &mut self.processors {
            proc.init();
        }
        self.initialized = true;
    }

    /// Pre-process all processors.
    pub fn pre_process_all(&mut self) {
        for proc in &mut self.processors {
            proc.pre_process();
        }
    }

    /// Process all processors that need processing at this timestep.
    pub fn process_all(
        &mut self,
        e_field: &VectorField3D,
        h_field: &VectorField3D,
        timestep: u64,
        time: f64,
    ) {
        for proc in &mut self.processors {
            if proc.needs_processing(timestep) {
                proc.process(e_field, h_field, timestep, time);
            }
        }
    }

    /// Post-process all processors.
    pub fn post_process_all(&mut self) {
        for proc in &mut self.processors {
            proc.post_process();
        }
    }

    /// Flush all processors.
    pub fn flush_all(&mut self) -> std::io::Result<()> {
        for proc in &mut self.processors {
            proc.flush()?;
        }
        Ok(())
    }

    /// Reset all processors.
    pub fn reset_all(&mut self) {
        for proc in &mut self.processors {
            proc.reset();
        }
        self.initialized = false;
    }

    /// Get a processor by index.
    pub fn get(&self, index: usize) -> Option<&dyn Processor> {
        self.processors.get(index).map(|p| p.as_ref())
    }

    /// Get a mutable reference to a boxed processor by index.
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Box<dyn Processor>> {
        self.processors.get_mut(index)
    }

    /// Get a processor by name.
    pub fn get_by_name(&self, name: &str) -> Option<&dyn Processor> {
        self.processors
            .iter()
            .find(|p| p.name() == name)
            .map(|p| p.as_ref())
    }

    /// Get a mutable reference to a boxed processor by name.
    pub fn get_by_name_mut(&mut self, name: &str) -> Option<&mut Box<dyn Processor>> {
        self.processors.iter_mut().find(|p| p.name() == name)
    }

    /// Remove a processor by index.
    pub fn remove(&mut self, index: usize) -> Option<Box<dyn Processor>> {
        if index < self.processors.len() {
            Some(self.processors.remove(index))
        } else {
            None
        }
    }

    /// Remove all processors.
    pub fn clear(&mut self) {
        self.processors.clear();
        self.initialized = false;
    }

    /// Iterate over processors.
    pub fn iter(&self) -> impl Iterator<Item = &dyn Processor> + '_ {
        self.processors.iter().map(|p| p.as_ref())
    }

    /// Iterate over boxed processors mutably.
    /// Use this to access processors with mutable references.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Box<dyn Processor>> {
        self.processors.iter_mut()
    }

    /// Get processor names.
    pub fn names(&self) -> Vec<&str> {
        self.processors.iter().map(|p| p.name()).collect()
    }

    /// Downcast and get a specific processor type.
    pub fn get_typed<T: Processor + 'static>(&self, index: usize) -> Option<&T> {
        self.processors
            .get(index)
            .and_then(|p| p.as_any().downcast_ref::<T>())
    }

    /// Downcast and get a specific processor type mutably.
    pub fn get_typed_mut<T: Processor + 'static>(&mut self, index: usize) -> Option<&mut T> {
        self.processors
            .get_mut(index)
            .and_then(|p| p.as_any_mut().downcast_mut::<T>())
    }
}

impl Default for ProcessingArray {
    fn default() -> Self {
        Self::new()
    }
}

/// Simple example processor for recording field values.
#[derive(Debug)]
pub struct FieldRecorder {
    name: String,
    position: [usize; 3],
    e_values: Vec<[f64; 3]>,
    h_values: Vec<[f64; 3]>,
    times: Vec<f64>,
    interval: u64,
}

impl FieldRecorder {
    /// Create a new field recorder.
    pub fn new(name: &str, position: [usize; 3]) -> Self {
        Self {
            name: name.to_string(),
            position,
            e_values: Vec::new(),
            h_values: Vec::new(),
            times: Vec::new(),
            interval: 1,
        }
    }

    /// Set recording interval.
    pub fn with_interval(mut self, interval: u64) -> Self {
        self.interval = interval.max(1);
        self
    }

    /// Get recorded E-field values.
    pub fn e_values(&self) -> &[[f64; 3]] {
        &self.e_values
    }

    /// Get recorded H-field values.
    pub fn h_values(&self) -> &[[f64; 3]] {
        &self.h_values
    }

    /// Get recorded times.
    pub fn times(&self) -> &[f64] {
        &self.times
    }
}

impl Processor for FieldRecorder {
    fn name(&self) -> &str {
        &self.name
    }

    fn process(
        &mut self,
        e_field: &VectorField3D,
        h_field: &VectorField3D,
        _timestep: u64,
        time: f64,
    ) {
        let [i, j, k] = self.position;

        self.e_values.push([
            e_field.x.get(i, j, k) as f64,
            e_field.y.get(i, j, k) as f64,
            e_field.z.get(i, j, k) as f64,
        ]);

        self.h_values.push([
            h_field.x.get(i, j, k) as f64,
            h_field.y.get(i, j, k) as f64,
            h_field.z.get(i, j, k) as f64,
        ]);

        self.times.push(time);
    }

    fn needs_processing(&self, timestep: u64) -> bool {
        self.interval == 0 || timestep.is_multiple_of(self.interval)
    }

    fn interval(&self) -> u64 {
        self.interval
    }

    fn reset(&mut self) {
        self.e_values.clear();
        self.h_values.clear();
        self.times.clear();
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn as_any_mut(&mut self) -> &mut dyn Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::arrays::Dimensions;

    #[test]
    fn test_processing_array_creation() {
        let array = ProcessingArray::new();
        assert!(array.is_empty());
    }

    #[test]
    fn test_add_processor() {
        let mut array = ProcessingArray::new();
        let recorder = FieldRecorder::new("test", [5, 5, 5]);

        array.add(recorder);
        assert_eq!(array.len(), 1);
    }

    #[test]
    fn test_process_all() {
        let mut array = ProcessingArray::new();
        array.add(FieldRecorder::new("recorder1", [2, 2, 2]));
        array.add(FieldRecorder::new("recorder2", [3, 3, 3]));

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        e_field.x.set(2, 2, 2, 1.0);
        e_field.x.set(3, 3, 3, 2.0);

        array.init_all();
        array.process_all(&e_field, &h_field, 0, 0.0);
        array.process_all(&e_field, &h_field, 1, 1e-12);

        // Check that values were recorded
        let recorder1 = array.get_typed::<FieldRecorder>(0).unwrap();
        assert_eq!(recorder1.times().len(), 2);
    }

    #[test]
    fn test_get_by_name() {
        let mut array = ProcessingArray::new();
        array.add(FieldRecorder::new("test_recorder", [0, 0, 0]));

        let proc = array.get_by_name("test_recorder");
        assert!(proc.is_some());
        assert_eq!(proc.unwrap().name(), "test_recorder");

        assert!(array.get_by_name("nonexistent").is_none());
    }

    #[test]
    fn test_field_recorder() {
        let mut recorder = FieldRecorder::new("test", [5, 5, 5]);

        let dims = Dimensions::new(10, 10, 10);
        let mut e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        e_field.x.set(5, 5, 5, 1.0);
        e_field.y.set(5, 5, 5, 2.0);
        e_field.z.set(5, 5, 5, 3.0);

        recorder.process(&e_field, &h_field, 0, 0.0);

        assert_eq!(recorder.e_values().len(), 1);
        assert!((recorder.e_values()[0][0] - 1.0).abs() < 1e-6);
        assert!((recorder.e_values()[0][1] - 2.0).abs() < 1e-6);
        assert!((recorder.e_values()[0][2] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_interval_processing() {
        let recorder = FieldRecorder::new("test", [0, 0, 0]).with_interval(10);

        assert!(recorder.needs_processing(0));
        assert!(!recorder.needs_processing(5));
        assert!(recorder.needs_processing(10));
        assert!(recorder.needs_processing(100));
    }

    #[test]
    fn test_reset() {
        let mut array = ProcessingArray::new();
        array.add(FieldRecorder::new("test", [0, 0, 0]));

        let dims = Dimensions::new(10, 10, 10);
        let e_field = VectorField3D::new(dims);
        let h_field = VectorField3D::new(dims);

        array.init_all();
        array.process_all(&e_field, &h_field, 0, 0.0);
        array.reset_all();

        let recorder = array.get_typed::<FieldRecorder>(0).unwrap();
        assert!(recorder.times().is_empty());
    }

    #[test]
    fn test_names() {
        let mut array = ProcessingArray::new();
        array.add(FieldRecorder::new("recorder1", [0, 0, 0]));
        array.add(FieldRecorder::new("recorder2", [1, 1, 1]));

        let names = array.names();
        assert_eq!(names.len(), 2);
        assert!(names.contains(&"recorder1"));
        assert!(names.contains(&"recorder2"));
    }
}
