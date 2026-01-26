use crate::arrays::{Dimensions, VectorField3D};
use crate::fdtd::operator::Operator;
use std::borrow::Cow;
use wgpu::util::DeviceExt;

/// GPU-accelerated FDTD engine using WebGPU.
pub struct GpuEngine {
    instance: wgpu::Instance,
    device: wgpu::Device,
    queue: wgpu::Queue,
    update_h_pipeline: wgpu::ComputePipeline,
    update_e_pipeline: wgpu::ComputePipeline,
    bind_group: wgpu::BindGroup,

    // Buffers
    e_field_buffer: wgpu::Buffer,
    h_field_buffer: wgpu::Buffer,

    // Dimensions
    dims: Dimensions,

    // Workgroup dispatch size
    dispatch_x: u32,
    dispatch_y: u32,
    dispatch_z: u32,
}

impl GpuEngine {
    /// Create a new GPU engine.
    pub fn new(operator: &Operator) -> Self {
        let dims = operator.dimensions();
        let total = dims.total();

        let instance = wgpu::Instance::default();
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
        }))
        .expect("Failed to find an appropriate adapter");

        let limits = adapter.limits();
        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: None,
            required_features: wgpu::Features::empty(),
            required_limits: limits,
            memory_hints: wgpu::MemoryHints::Performance,
            ..Default::default()
        }))
        .expect("Failed to create device");

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("FDTD Shader"),
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(include_str!("shaders.wgsl"))),
        });

        let field_size = (total * 3 * std::mem::size_of::<f32>()) as u64;

        let e_field_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("E Field Buffer"),
            size: field_size,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let h_field_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("H Field Buffer"),
            size: field_size,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_DST
                | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let e_coeff = operator.e_coefficients();
        let h_coeff = operator.h_coefficients();

        let mut e_coeff_data: Vec<f32> = Vec::with_capacity(total * 6);
        for i in 0..3 {
            e_coeff_data.extend_from_slice(e_coeff.ca[i].as_slice());
        }
        for i in 0..3 {
            e_coeff_data.extend_from_slice(e_coeff.cb[i].as_slice());
        }

        let mut h_coeff_data: Vec<f32> = Vec::with_capacity(total * 6);
        for i in 0..3 {
            h_coeff_data.extend_from_slice(h_coeff.da[i].as_slice());
        }
        for i in 0..3 {
            h_coeff_data.extend_from_slice(h_coeff.db[i].as_slice());
        }

        let e_coeff_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("E Coeff Buffer"),
            contents: bytemuck::cast_slice(&e_coeff_data),
            usage: wgpu::BufferUsages::STORAGE,
        });

        let h_coeff_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("H Coeff Buffer"),
            contents: bytemuck::cast_slice(&h_coeff_data),
            usage: wgpu::BufferUsages::STORAGE,
        });

        // Uniform buffer removed, constants used instead

        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("FDTD Bind Group Layout"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: false },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 3,
                    visibility: wgpu::ShaderStages::COMPUTE,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Storage { read_only: true },
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
            ],
        });

        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("FDTD Bind Group"),
            layout: &bind_group_layout,
            entries: &[
                wgpu::BindGroupEntry {
                    binding: 0,
                    resource: e_field_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: h_field_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 2,
                    resource: e_coeff_buffer.as_entire_binding(),
                },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: h_coeff_buffer.as_entire_binding(),
                },
            ],
        });

        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("Pipeline Layout"),
            bind_group_layouts: &[&bind_group_layout],
            ..Default::default()
        });

        let constants_data = [
            ("nx", dims.nx as f64),
            ("ny", dims.ny as f64),
            ("nz", dims.nz as f64),
        ];

        let update_h_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Update H Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("update_h"),
            compilation_options: wgpu::PipelineCompilationOptions {
                constants: &constants_data,
                ..Default::default()
            },
            cache: None,
        });

        let update_e_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some("Update E Pipeline"),
            layout: Some(&pipeline_layout),
            module: &shader,
            entry_point: Some("update_e"),
            compilation_options: wgpu::PipelineCompilationOptions {
                constants: &constants_data,
                ..Default::default()
            },
            cache: None,
        });

        // Dispatch dimensions swapped to match shader mapping:
        // GlobalID.x -> k (fastest mem) -> dispatch_x covers nz
        // GlobalID.y -> j -> dispatch_y covers ny
        // GlobalID.z -> i -> dispatch_z covers nx
        let dispatch_x = (dims.nz as u32 + 63) / 64; // Workgroup size 64 along X
        let dispatch_y = (dims.ny as u32 + 1) / 2; // Workgroup size 2 along Y
        let dispatch_z = (dims.nx as u32 + 1) / 2; // Workgroup size 2 along Z

        Self {
            instance,
            device,
            queue,
            update_h_pipeline,
            update_e_pipeline,
            bind_group,
            e_field_buffer,
            h_field_buffer,
            dims,
            dispatch_x,
            dispatch_y,
            dispatch_z,
        }
    }

    /// Perform one FDTD timestep on the GPU.
    pub fn step(&self) {
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Step Encoder"),
            });

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Update H Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.update_h_pipeline);
            compute_pass.set_bind_group(0, &self.bind_group, &[]);
            compute_pass.dispatch_workgroups(self.dispatch_x, self.dispatch_y, self.dispatch_z);
        }

        {
            let mut compute_pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: Some("Update E Pass"),
                timestamp_writes: None,
            });
            compute_pass.set_pipeline(&self.update_e_pipeline);
            compute_pass.set_bind_group(0, &self.bind_group, &[]);
            compute_pass.dispatch_workgroups(self.dispatch_x, self.dispatch_y, self.dispatch_z);
        }

        self.queue.submit(Some(encoder.finish()));
    }

    /// Upload E-field data to the GPU.
    pub fn write_e_field(&self, field: &VectorField3D) {
        let mut data = Vec::with_capacity(self.dims.total() * 3);
        data.extend_from_slice(field.x.as_slice());
        data.extend_from_slice(field.y.as_slice());
        data.extend_from_slice(field.z.as_slice());

        self.queue
            .write_buffer(&self.e_field_buffer, 0, bytemuck::cast_slice(&data));
    }

    /// Upload H-field data to the GPU.
    pub fn write_h_field(&self, field: &VectorField3D) {
        let mut data = Vec::with_capacity(self.dims.total() * 3);
        data.extend_from_slice(field.x.as_slice());
        data.extend_from_slice(field.y.as_slice());
        data.extend_from_slice(field.z.as_slice());

        self.queue
            .write_buffer(&self.h_field_buffer, 0, bytemuck::cast_slice(&data));
    }

    /// Download E-field data from the GPU.
    pub fn read_e_field(&self, field: &mut VectorField3D) {
        self.read_buffer_to_field(&self.e_field_buffer, field);
    }

    /// Download H-field data from the GPU.
    pub fn read_h_field(&self, field: &mut VectorField3D) {
        self.read_buffer_to_field(&self.h_field_buffer, field);
    }

    /// Update a single element of the E-field on the GPU.
    pub fn update_e_field_element(&self, i: usize, j: usize, k: usize, dir: usize, value: f32) {
        let total = self.dims.total();
        let offset_elems = match dir {
            0 => 0,
            1 => total,
            2 => 2 * total,
            _ => return,
        };
        let idx = self.dims.to_linear(i, j, k);
        let final_offset = (offset_elems + idx) * 4;

        self.queue.write_buffer(
            &self.e_field_buffer,
            final_offset as u64,
            bytemuck::bytes_of(&value),
        );
    }

    fn read_buffer_to_field(&self, buffer: &wgpu::Buffer, field: &mut VectorField3D) {
        let size = buffer.size();
        let staging_buffer = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("Staging Buffer"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor {
                label: Some("Read Encoder"),
            });

        encoder.copy_buffer_to_buffer(buffer, 0, &staging_buffer, 0, size);
        self.queue.submit(Some(encoder.finish()));

        let slice = staging_buffer.slice(..);
        let (sender, receiver) = futures::channel::oneshot::channel();
        slice.map_async(wgpu::MapMode::Read, move |v| sender.send(v).unwrap());

        self.instance.poll_all(true);
        pollster::block_on(receiver).unwrap().unwrap();

        let data = slice.get_mapped_range();
        let floats: &[f32] = bytemuck::cast_slice(&data);

        let total = self.dims.total();
        field.x.as_mut_slice().copy_from_slice(&floats[0..total]);
        field
            .y
            .as_mut_slice()
            .copy_from_slice(&floats[total..2 * total]);
        field
            .z
            .as_mut_slice()
            .copy_from_slice(&floats[2 * total..3 * total]);

        drop(data);
        staging_buffer.unmap();
    }

    /// Reset all fields to zero.
    pub fn reset(&self) {
        let size = self.e_field_buffer.size();
        let mut encoder = self
            .device
            .create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.clear_buffer(&self.e_field_buffer, 0, Some(size));
        encoder.clear_buffer(&self.h_field_buffer, 0, Some(size));
        self.queue.submit(Some(encoder.finish()));
    }

    /// Wait for all GPU operations to complete.
    pub fn wait_idle(&self) {
        self.instance.poll_all(true);
    }
}
