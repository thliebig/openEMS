%
% Tutorials / 7T MRI Loop Coil
%
% Tested with
%  - Octave 11.3
%  - openEMS v0.37
%
% (C) 2013-2026 Thorsten Liebig <thorsten.liebig@gmx.de>

close all
clear
clc

%% Simulation Parameters
%% ---------------------
%% This tutorial models a surface loop coil designed for 7 T MRI (proton
%% Larmor frequency 298 MHz) placed next to the "Ella" Virtual Family voxel
%% body model. The loop is tuned to resonance with lumped capacitors;
%% ``loop.C_gap`` and ``loop.port_R`` set the tuning capacitance and the
%% feed resistance that represents the transmit/receive switch. All
%% dimensions are in millimetres.
physical_constants; %get some physical constants like c0 and MUE0
unit = 1e-3; % all length in mm

% Loop-Coil parameter
loop.length = 80;       % length of the loop (in z-direction)
loop.width = 60;        % width of the loop  (in y-direction)
loop.strip_width = 5;   % metal strip width
loop.strip_N_cells = 3; % number of cells over the strip length
loop.air_gap = loop.strip_width/3;       % air gap width for lumped capacitors
loop.pos_x = -130;       % position of loop
loop.C_gap = 5.4e-12;   % lumped cap value
loop.port_R = 2.5;      % feeding port resistance

%% Human Body Model Setup
%% ----------------------
%% The Virtual Family voxel model provides realistic, frequency-dependent
%% dielectric tissue properties derived from measured data. ``Convert_VF_DiscMaterial``
%% slices the full-body dataset to the head/shoulder region (``body_model_range``)
%% at 298 MHz, which keeps memory use manageable. This conversion only needs
%% to run once; subsequent script runs skip it if the HDF5 output file
%% already exists. If the VF dataset is not installed the script falls back
%% automatically to a simple ellipsoidal phantom with average head-tissue
%% properties at 298 MHz — no manual changes are required.
% set file name for human body model to create with "Convert_VF_DiscMaterial"
% the file name should contain a full path
body_model_file = [pwd '/Ella_centered_298MHz.h5'];

% convert only part of the model (head/shoulder section)
body_model_range = {[],[],[-0.85 -0.4]};

% paths to virtual family voxel models (VFVM), adept to your install!
VF_raw_filesuffix = '/tmp/Ella_26y_V2_1mm';
VF_mat_db_file = '/tmp/DB_h5_20120711_SEMCADv14.8.h5';

% delete(body_model_file); % uncomment to delete old model if something changed

% use cached HDF5 if it already exists; otherwise try to convert from raw VF files;
% fall back to a homogeneous phantom if the VF dataset is unavailable
use_body_model = exist(body_model_file, 'file') == 2;
if ~use_body_model
    try
        Convert_VF_DiscMaterial(VF_raw_filesuffix, VF_mat_db_file, body_model_file, ...
                                'Frequency', 298e6, 'Center', 1, ...
                                'Range', body_model_range);
        use_body_model = 1;
    catch
        warning('openEMS:MRI_Loop_Coil', ...
            'VF body model not found — using homogeneous ellipsoidal phantom fallback.');
    end
end

% rotate model to face the nose in x-dir, and translate
body_model_transform = {'Rotate_X',pi,'Rotate_Z',pi/2, ...
                        'Translate',[0,5,-720]};

% the head should + part of shoulder should fit this box
body_box.start = [-120 -150 -200];
body_box.stop  = [+100 +150 +130];

% box with high res mesh
mesh_box.start = [-120 -80 -120];
mesh_box.stop  = [+100 +80 +120];
mesh_box.resolution = 2;

%% Air Box Size
%% ------------
%% The air box extends the simulation domain beyond the structure so that
%% the absorbing boundary conditions operate in free space rather than
%% inside heterogeneous tissue. A 150 mm margin is sufficient at 298 MHz
%% to keep reflections from the MUR boundaries below the solver's
%% convergence threshold before they can reach the coil.
Air_Box = 150;      % size of the surrounding air box (150mm)

%% FDTD Solver and Excitation
%% --------------------------
%% A Gaussian pulse centred at 298 MHz with a 300 MHz corner frequency
%% excites all resonances in a single time-domain run, from which
%% frequency-domain results are extracted by DFT. ``CellConstantMaterial``
%% disables the default sub-cell material averaging so that every Yee cell
%% holds a single homogeneous material. That models the tissue boundaries
%% more coarsely, but it is what the SAR calculation assumes and what
%% IEC/IEEE 62704-1 requires. MUR first-order absorbing boundaries
%% terminate the domain on all six faces.
% init FDTD structure
FDTD = InitFDTD( 'EndCriteria', 1e-4, 'CellConstantMaterial', 1); % constant material per voxel, required for SAR

% define gaussian pulse excitation signal
f0 = 298e6; % center frequency
fc = 300e6; % 20 dB corner frequency
FDTD = SetGaussExcite( FDTD, f0, fc );

% setup boundary conditions
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% Geometry Initialisation
%% -----------------------
%% CSXCAD holds all geometric and material objects that openEMS translates
%% into FDTD coefficients. Initialising the empty structure here provides
%% the container into which loop metal, lumped elements, the body model,
%% and field dump boxes are added in the sections below.
CSX = InitCSX();

%% Loop Coil Geometry
%% ------------------
%% The loop consists of flat copper strip segments interrupted at four gaps
%% where lumped capacitors provide distributed resonance tuning. Three gaps
%% carry tuning capacitors (``caps_y`` and ``caps_z``); the fourth gap holds
%% the lumped feed port whose series resistance ``loop.port_R`` models the
%% transmit/receive switch or preamplifier input. The strip width and gap
%% size together set the capacitor footprint and self-resonance frequency.
% setup all properties needed
CSX = AddMetal( CSX, 'loop' );
CSX = AddLumpedElement( CSX, 'caps_y', 1, 'C', loop.C_gap);
CSX = AddLumpedElement( CSX, 'caps_z', 2, 'C', loop.C_gap);

% horizontal (y-direction) strips
start = [loop.pos_x -loop.width/2   -loop.length/2];
stop  = [loop.pos_x -loop.air_gap/2 -loop.length/2+loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x  -loop.width/2   loop.length/2 ];
stop  = [loop.pos_x  -loop.air_gap/2 loop.length/2-loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2   -loop.length/2];
stop  = [loop.pos_x loop.air_gap/2 -loop.length/2+loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2   loop.length/2  ];
stop  = [loop.pos_x loop.air_gap/2 loop.length/2-loop.strip_width];
CSX = AddBox(CSX,'loop',10,start,stop);

% vertical (z-direction) strips
start = [loop.pos_x -loop.width/2                  -loop.length/2+loop.strip_width];
stop  = [loop.pos_x -loop.width/2+loop.strip_width -loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x -loop.width/2                  loop.length/2-loop.strip_width];
stop  = [loop.pos_x -loop.width/2+loop.strip_width loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2                   -loop.length/2+loop.strip_width];
stop  = [loop.pos_x loop.width/2-loop.strip_width  -loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

start = [loop.pos_x loop.width/2                   loop.length/2-loop.strip_width ];
stop  = [loop.pos_x loop.width/2-loop.strip_width  loop.air_gap/2];
CSX = AddBox(CSX,'loop',10,start,stop);

% add the lumped capacities
start = [loop.pos_x -loop.width/2+loop.strip_width/2-loop.air_gap/2 -loop.air_gap/2];
stop  = [loop.pos_x -loop.width/2+loop.strip_width/2+loop.air_gap/2 +loop.air_gap/2];
CSX = AddBox(CSX,'caps_z',10,start,stop);

start = [loop.pos_x loop.width/2-loop.strip_width/2-loop.air_gap/2 -loop.air_gap/2];
stop  = [loop.pos_x loop.width/2-loop.strip_width/2+loop.air_gap/2 +loop.air_gap/2];
CSX = AddBox(CSX,'caps_z',10,start,stop);

start = [loop.pos_x -loop.air_gap/2 loop.length/2-loop.strip_width/2-loop.air_gap/2];
stop  = [loop.pos_x +loop.air_gap/2 loop.length/2-loop.strip_width/2+loop.air_gap/2];
CSX = AddBox(CSX,'caps_y',10,start,stop);

% add a lumped port as excitation
start = [loop.pos_x -loop.air_gap/2 -loop.length/2+loop.strip_width/2-loop.air_gap/2];
stop  = [loop.pos_x +loop.air_gap/2 -loop.length/2+loop.strip_width/2+loop.air_gap/2];
[CSX port] = AddLumpedPort(CSX, 100, 1, loop.port_R, start, stop, [0 1 0], true);

%% Body Model Placement
%% --------------------
%% When the VF dataset is available, the pre-converted HDF5 discrete material
%% file is linked into CSXCAD and bounded by ``body_box``, which crops the
%% dataset to the shoulder-to-crown region; ``Scale`` and ``Transform``
%% orient the voxel grid so that Ella faces in the x-direction with the head
%% centre positioned correctly relative to the loop at x = -130 mm. Without
%% VF data a bundled 3-layer ellipsoidal head phantom (skin / skull / brain,
%% tissue properties at 298 MHz from the IT'IS database) is used instead.
%% It produces the same SAR and B1 workflow output with a realistic three-tissue
%% distribution; the same ``Scale`` applies as the phantom mesh is in metres.
if use_body_model
    CSX = AddDiscMaterial(CSX, 'body_model', 'File', body_model_file, 'Scale', 1/unit, 'Transform', body_model_transform);
    CSX = AddBox(CSX, 'body_model', 0, body_box.start, body_box.stop);
else
    phantom_file = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'resources', 'phantoms', 'phantom_head_298MHz.h5');
    CSX = AddDiscMaterial(CSX, 'body_model', 'File', phantom_file, 'Scale', 1/unit);
    CSX = AddBox(CSX, 'body_model', 0, body_box.start, body_box.stop);
end

%% Mesh Generation
%% ---------------
%% ``DetectEdges`` seeds mesh lines at every conductor boundary and material
%% interface; these are refined to 2 mm inside the body region
%% (``mesh_box.resolution``) to resolve tissue heterogeneity accurately.
%% The outer domain is graded more coarsely — up to 10 cells per wavelength
%% — to limit the total cell count, and the 150 mm air spacers give the
%% absorbing boundaries enough room to attenuate outgoing waves.
% create loop mesh
mesh = DetectEdges(CSX);

% add a dense homogeneous mesh inside the human body model
mesh.x = [mesh.x mesh_box.start(1) mesh_box.stop(1)];
mesh.y = [mesh.y mesh_box.start(2) mesh_box.stop(2)];
mesh.z = [mesh.z mesh_box.start(3) mesh_box.stop(3)];

% add lines in x-dir for the loop and a cell centered around 0
mesh.x = [mesh.x loop.pos_x -mesh_box.resolution/2 mesh_box.resolution/2];

% smooth the mesh for the loop & body
mesh = SmoothMesh(mesh, mesh_box.resolution);

% add air spacer
mesh.x = [-Air_Box+mesh.x(1) mesh.x mesh.x(end)+Air_Box];
mesh.y = [-Air_Box+mesh.y(1) mesh.y mesh.y(end)+Air_Box];
mesh.z = [-Air_Box+mesh.z(1) mesh.z mesh.z(end)+Air_Box];

mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 40, 1.5, 'algorithm', 1);

%% Field and SAR Dump Boxes
%% ------------------------
%% Frequency-domain H-field (DumpType 11) and local SAR (DumpType 20,
%% DumpMode 2) are recorded on the axial (xy) and sagittal (xz) planes
%% through the head centre. Thin 2D slices keep output file sizes small
%% while capturing the clinically important field distributions needed to
%% assess transmit homogeneity and specific absorption rate compliance.
CSX = AddDump(CSX,'Hf_xy','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xy',0, body_box.start.*[1 1 0], body_box.stop.*[1 1 0]);
CSX = AddDump(CSX,'SAR_xy','DumpType',20,'DumpMode',2,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'SAR_xy',0, body_box.start.*[1 1 0], body_box.stop.*[1 1 0]);

CSX = AddDump(CSX,'Hf_xz','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xz',0, body_box.start.*[1 0 1], body_box.stop.*[1 0 1]);
CSX = AddDump(CSX,'SAR_xz','DumpType',20,'DumpMode',2,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'SAR_xz',0, body_box.start.*[1 0 1], body_box.stop.*[1 0 1]);

%% Absorbing Boundary Padding
%% --------------------------
%% ``AddPML`` inserts ten extra mesh lines at each face of the simulation
%% domain to host the perfectly matched layer or MUR absorbing boundaries
%% without distorting the interior mesh. These lines are added after mesh
%% smoothing so that they do not interfere with the automatic grading
%% applied inside the body region.
mesh = AddPML(mesh, 10);

%% Finalise and Apply Mesh
%% -----------------------
%% ``DefineRectGrid`` commits the completed mesh line vectors to the CSXCAD
%% structure, converting them to the Yee-cell grid that openEMS uses
%% internally. Printing the cell count beforehand gives a quick estimate of
%% expected run time and memory demand before launching the solver.
disp(['number of cells: ' num2str(1e-6*numel(mesh.x)*numel(mesh.y)*numel(mesh.z)) ' Mcells'])
CSX = DefineRectGrid( CSX, unit, mesh );

%% Simulation Folder Preparation
%% ------------------------------
%% All openEMS input and output files are written to a dedicated
%% subdirectory named after the script. ``CleanupSimPath`` removes any
%% previous results so that stale data from an earlier run cannot
%% contaminate the current post-processing.
Sim_Path = ['tmp_' mfilename];
Sim_CSX = [mfilename '.xml'];

CleanupSimPath(Sim_Path);

%% Write Simulation XML
%% --------------------
%% ``WriteOpenEMS`` serialises the FDTD settings and the CSXCAD geometry
%% to an XML file, which is the sole input openEMS reads at startup.
%% Keeping this write step separate from the solver call makes it easy
%% to inspect the setup or re-run the solver without re-executing the
%% full geometry script.
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% Visualise and Export Geometry
%% ------------------------------
%% ``CSXGeomPlot`` renders the geometry for a visual sanity check and
%% simultaneously exports VTK polydata for Paraview. The
%% ``--RenderDiscMaterial`` flag makes the body voxel model visible
%% alongside the coil metal, which is essential for confirming that the
%% loop is positioned correctly relative to the tissue before committing
%% to a long solver run.
CSXGeomPlot( [Sim_Path '/' Sim_CSX] , ['--export-polydata-vtk=' Sim_Path ' --RenderDiscMaterial -v']);

%% Run the Simulation
%% ------------------
%% ``RunOpenEMS`` launches the FDTD solver, which reads the XML file and
%% advances in time until the stored energy in the domain falls below the
%% ``EndCriteria`` threshold (-40 dB relative to peak). On a typical
%% workstation this problem takes several hours owing to the fine mesh
%% required to resolve heterogeneous tissue at 298 MHz.
RunOpenEMS( Sim_Path, Sim_CSX);

%% Port Post-Processing
%% --------------------
%% ``calcPort`` computes the complex voltage and current spectra at the
%% feed port by DFT of the time-domain recordings. The input impedance
%% ``Zin``, the reflection coefficient ``s11``, and the accepted power
%% ``P0_in`` at 298 MHz are derived here; ``P0_in`` serves as the
%% normalisation factor for all subsequent field and SAR quantities.
freq = linspace( f0-fc, f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

% get the feeding power for frequency f0
P0_in = interp1(freq, port.P_acc, f0);

%% S-Parameter and Admittance Plots
%% ---------------------------------
%% S11 shows whether the coil is well matched at 298 MHz — a deep notch
%% at the target frequency indicates resonance with low reflection. The
%% input admittance plot separates resistive loss from reactive stored
%% energy and confirms that the capacitor values have cancelled the
%% reactive part at resonance.
% plot reflection coefficient S11
figure
h = plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}| (dB)' );

% plot feed point admittance
figure
h = plot( freq/1e6, real(1./Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(1./Zin), 'r--', 'Linewidth', 2 );
title( 'feed port admittance' );
xlabel( 'frequency f (MHz)' );
ylabel( 'admittance Y_{in} (S)' );
legend( 'real', 'imag' );

%% SAR Distribution -- Axial Plane (XY)
%% -------------------------------------
%% The absorbed power per unit mass (SAR) is normalised to 1 W accepted
%% input power so that results scale directly to any desired transmit
%% level. The axial slice at z = 0 reveals the hot-spot pattern caused
%% by wave interference inside the head, which is the primary safety
%% figure of merit for MRI coil assessment under IEC 60601-2-33.
[SAR SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR_xy.h5']);
SAR = SAR.FD.values{1}/P0_in;

% SAR plot
figure()
subplot(1,2,1);
[X Y] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
colormap('hot');
h = pcolor(X,Y,(squeeze(SAR)));
% h = pcolor(X,Y,log10(squeeze(SAR)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('local SAR');
axis equal tight

%% SAR Distribution -- Sagittal Plane (XZ)
%% -----------------------------------------
%% The sagittal (xz) cut complements the axial view by showing how SAR
%% varies with depth along the head-to-foot axis. Together the two planes
%% provide enough spatial context to identify the peak local SAR location,
%% which determines IEC safety compliance at the chosen transmit power.
[SAR SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR_xz.h5']);
SAR = SAR.FD.values{1}/P0_in;

% SAR plot
subplot(1,2,2);
[X Z] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{3});
colormap('hot');
h = pcolor(X,Z,(squeeze(SAR)));
% h = pcolor(X,Y,log10(squeeze(SAR)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('z -->');
title('local SAR');
axis equal tight

%% B1 Field Maps -- Axial Plane (XY)
%% -----------------------------------
%% B1+ = (Bx + j*By)/2 is the circularly polarised transmit field that
%% flips nuclear spins; B1- = (Bx - j*By)/2 is the receive sensitivity.
%% Both are normalised to the square root of accepted power (T/sqrt(W))
%% so they characterise the coil independently of drive level. The
%% logarithmic colour scale reveals the large dynamic range typical of
%% surface coils placed close to tissue.
[H_field H_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xy.h5']);
% calc Bx,By, B1p, B1m normalize to the input-power
Bx = MUE0*H_field.FD.values{1}(:,:,:,1)/sqrt(P0_in);
By = MUE0*H_field.FD.values{1}(:,:,:,2)/sqrt(P0_in);
B1p = 0.5*(Bx+1j*By);
B1m = 0.5*(Bx-1j*By);
% create a 2D grid to plot on
[X Y] = ndgrid(H_mesh.lines{1},H_mesh.lines{2});

Dump2VTK([Sim_Path '/B1p_xy.vtk'], abs(B1p), H_mesh, 'B-Field');
Dump2VTK([Sim_Path '/B1m_xy.vtk'], abs(B1m), H_mesh, 'B-Field');

% B1+ plot
figure()
subplot(1,2,1);
h = pcolor(X,Y,log10(abs(B1p)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^+ field (dB)');
axis equal tight

% B1- plot
subplot(1,2,2);
h = pcolor(X,Y,log10(abs(B1m)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^- field (dB)');
axis equal tight

%% B1 Field Maps -- Sagittal Plane (XZ)
%% --------------------------------------
%% The sagittal B1+/B1- maps show the depth penetration of transmit and
%% receive sensitivity from the loop, which is placed laterally at
%% x = -130 mm. The rapid roll-off with distance is characteristic of
%% surface coils and limits their useful field of view compared with
%% volume resonators such as birdcage coils.
[H_field H_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xz.h5']);
% calc Bx,By, B1p, B1m normalize to the input-power
Bx = MUE0*H_field.FD.values{1}(:,:,:,1)/sqrt(P0_in);
By = MUE0*H_field.FD.values{1}(:,:,:,2)/sqrt(P0_in);
B1p = 0.5*(Bx+1j*By);
B1m = 0.5*(Bx-1j*By);
% create a 2D grid to plot on
[X Z] = ndgrid(H_mesh.lines{1},H_mesh.lines{3});

Dump2VTK([Sim_Path '/B1p_xz.vtk'], abs(B1p), H_mesh, 'B-Field');
Dump2VTK([Sim_Path '/B1m_xz.vtk'], abs(B1m), H_mesh, 'B-Field');

% B1+ plot
figure()
subplot(1,2,1);
h = pcolor(X,Z,log10(squeeze(abs(B1p))));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('z -->');
title('B_1^+ field (dB)');
axis equal tight

% B1- plot
subplot(1,2,2);
h = pcolor(X,Z,log10(squeeze(abs(B1m))));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('z -->');
title('B_1^- field (dB)');
axis equal tight

%% Export SAR Data to VTK
%% -----------------------
%% ``ConvertHDF5_VTK`` rescales the raw HDF5 field dump by 1/P0_in and
%% writes VTK files that Paraview can overlay on the 3D body model. This
%% enables interactive exploration of the SAR distribution in three
%% dimensions, which is more informative than the 2D slice plots above.
ConvertHDF5_VTK([Sim_Path '/SAR_xy.h5'],[Sim_Path '/SAR_xy'], 'weight', 1/P0_in, 'FieldName', 'SAR');
ConvertHDF5_VTK([Sim_Path '/SAR_xz.h5'],[Sim_Path '/SAR_xz'], 'weight', 1/P0_in, 'FieldName', 'SAR');

%% Export B1 Field Data to VTK
%% ----------------------------
%% The H-field HDF5 dumps are converted to VTK format with a weighting
%% of MUE0/sqrt(P0_in) to yield the B1 field in T/sqrt(W). Loading
%% these files alongside the SAR data in Paraview allows the correlation
%% between high-field regions and elevated SAR to be visualised directly
%% in three dimensions.
ConvertHDF5_VTK([Sim_Path '/Hf_xy.h5'],[Sim_Path '/B1_xy'], 'weight', MUE0/sqrt(P0_in), 'FieldName', 'B1-field');
ConvertHDF5_VTK([Sim_Path '/Hf_xz.h5'],[Sim_Path '/B1_xz'], 'weight', MUE0/sqrt(P0_in), 'FieldName', 'B1-field');
