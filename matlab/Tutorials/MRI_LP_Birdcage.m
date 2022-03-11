%
% Tutorials / 3T MRI Low Pass Birdcage coil
%
% Description at:
% http://openems.de/index.php/Tutorial:_MRI_LP_Birdcage
%
% Estimated time to run:    ~7h @ ~65MC/s
% Memory requirement (RAM): ~ 700MB
%
% Tested with
%  - openEMS v0.0.33
%  - Matlab 7.12.0 (R2011a)
%
% (C) 2013-2015 Thorsten Liebig <thorsten.liebig@gmx.de>


close all
clear
clc

% simulation setup
f0 = 128e6;
excite.f_0 = 75e6; % excite gaussian pulse center frequency
excite.f_c = 75e6;  % excite gaussian pulse cutoff frequency

postproc_only = 0;  % set to 1 to perform only post processing
GeomPlot = 1;       % set to 0 to skip geometry viewer

% bore setup
Bore.rad = 320;
Bore.length = 1600;

% birdcage setup
BC.N_rungs = 8;
BC.rad = 120;
BC.stripwidth = 10;
BC.portwidth = BC.stripwidth/2;
BC.portlength = BC.stripwidth/2;
BC.length = 250;
BC.cap = 2.6e-12;

% feed amplitude and phase at given rungs
BC.feed_pos = [1 3];
BC.feed_amp = [1 -1j];

%% define the human body model (virtual family)
% set file name for human body model to create with "Convert_VF_DiscMaterial"
% the file name should contain a full path
body_model_file = [pwd '/Ella_centered_' num2str(f0/1e6) 'MHz.h5'];

% convert only part of the model (head/shoulder section)
body_model_range = {[],[],[-0.85 0]};

body_mesh_res = 2.5; % should be something like: BC.stripwidth/4

% paths to virtual family voxel models (VFVM), adept to your install!
VF_raw_filesuffix = '/tmp/Ella_26y_V2_1mm';
VF_mat_db_file = '/tmp/DB_h5_20120711_SEMCADv14.8.h5';

% delete(body_model_file); % uncomment to delete old model if something changed

% convert model (if it does not exist)
Convert_VF_DiscMaterial(VF_raw_filesuffix, VF_mat_db_file, body_model_file, ...
                        'Frequency', f0, 'Center', 1, ...
                        'Range', body_model_range);

% rotate model to face the nose in +y-dir, and translate
body_model_transform = {'Rotate_X',pi,'Rotate_Z',pi, ...
                        'Translate',[0,5,-720]};

%% some internal parameter
physical_constants % load important physical constants
end_crit = 1e-5;    %abort simulation at -50dB energy drop
unit = 1e-3;        %drawing unit used

%capacity footprint is 4mm x 4mm
lambda_min = c0/(excite.f_0+excite.f_c);

% meshing options
% desired mesh resolution
mesh_res([1 3]) = min(15,lambda_min/20/unit);
mesh_res(2) = body_mesh_res / BC.rad;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD('CoordSystem', 1, ... %init a cylindrical FDTD setup
    'EndCriteria', 1e-4, ... % with an end criteria of -40dB (1e-4)
    'MultiGrid', '10,20',... % add two cylindrical sub-grids at a radius of 10 and 20 mm
    'CellConstantMaterial', 1); % assume a material is constant inside
                                % a cell (material probing in cell center)

% define the excitation time-signal (unmodulated gaussian pulse)
FDTD = SetGaussExcite(FDTD,excite.f_0,excite.f_c);

% define & set boundary conditions
%   - pml in +/- z-direction
%   - boundaries in -r and +/- alpha direction disabled (full cylindrical mesh)
%   - PEC boundary in +r-direction to model bore RF shield
FDTD = SetBoundaryCond(FDTD, [0 0 0 0 3 3]);


%% setup CSXCAD geometry & mesh (cylindrical)
CSX = InitCSX('CoordSystem',1);

% init empty mesh structure
mesh.r = [];
mesh.a = [];
mesh.z = [];

%% Create metal bird cage and rung capacities
CSX = AddMetal(CSX,'metal');
CSX = AddLumpedElement(CSX,'caps','z','C',BC.cap);

da_Strip = BC.stripwidth/BC.rad; % width of a strip in radiant
da_Caps = BC.portwidth/BC.rad;   % width of a cap/port in radiant
da_Segs = 2*pi/BC.N_rungs;       % width of a rung in radiant

a_start = -pi-da_Segs/2;         % starting angle

w0 = 2*pi*f0;
T0 = 1/f0;

% port counter
port_Nr = 1;

a0 = a_start;

for n=1:BC.N_rungs
    start = [BC.rad a0+da_Segs/2-da_Caps/2 -0.5*BC.portlength];
    stop  = [BC.rad a0+da_Segs/2+da_Caps/2 +0.5*BC.portlength];
    CSX = AddBox(CSX,'caps',1, start, stop);
    
    start = [BC.rad a0+da_Segs/2-da_Caps/2 0.5*BC.length-BC.stripwidth/2-BC.portlength];
    stop  = [BC.rad a0+da_Segs/2+da_Caps/2 0.5*BC.length-BC.stripwidth/2];
    if (~isempty(intersect(n, BC.feed_pos)) && (BC.feed_amp(port_Nr)~=0)) % active port
        exc_amp = abs(BC.feed_amp(port_Nr));
        
        % calculate time delay to achieve a given phase shift at f0
        T = -angle(BC.feed_amp(port_Nr)) / w0;
        if T<0
            T = T + T0;
        end
        [CSX port{port_Nr}] = AddLumpedPort(CSX, 100, port_Nr, 50, start, stop, [0 0 1]*exc_amp, true,'Delay',T);
        
        %increase port count
        port_Nr = port_Nr+1;
        
        start = [BC.rad a0+da_Segs/2-da_Strip/2 0.5*BC.length-BC.stripwidth/2-BC.portlength];
    elseif ~isempty(intersect(n, BC.feed_pos))  % passive port
        [CSX port{port_Nr}] = AddLumpedPort(CSX, 100, port_Nr, 50, start, stop, [0 0 1], false);
        
        %increase port count
        port_Nr = port_Nr+1;
        
        start = [BC.rad a0+da_Segs/2-da_Strip/2 0.5*BC.length-BC.stripwidth/2-BC.portlength];
    else
        start = [BC.rad a0+da_Segs/2-da_Strip/2 0.5*BC.length];
    end
    
    % the start z-coordinate depends on the port (see above)
    stop  = [BC.rad a0+da_Segs/2+da_Strip/2 0.5*BC.portlength];
    CSX = AddBox(CSX,'metal',1, start, stop);
    
    start = [BC.rad a0+da_Segs/2-da_Strip/2 -0.5*BC.length];
    stop  = [BC.rad a0+da_Segs/2+da_Strip/2 -0.5*BC.portlength];
    CSX = AddBox(CSX,'metal',1, start, stop);
    
    % some additional mesh lines
    mesh.a = [mesh.a a0+da_Segs/2];
    
    a0 = a0 + da_Segs;
end

% create metal top ring
start = [BC.rad a_start      -(BC.length-BC.stripwidth)/2];
stop  = [BC.rad a_start+2*pi -(BC.length+BC.stripwidth)/2];
CSX = AddBox(CSX,'metal',1, start, stop);

% create metal bottom ring
start = [BC.rad a_start      (BC.length-BC.stripwidth)/2];
stop  = [BC.rad a_start+2*pi (BC.length+BC.stripwidth)/2];
CSX = AddBox(CSX,'metal',1, start, stop);

%% create smooth mesh
mesh = DetectEdges(CSX, mesh);
mesh.r = [0 SmoothMeshLines([body_mesh_res*1.5 mesh.r], body_mesh_res)];
mesh.z = SmoothMeshLines(mesh.z, body_mesh_res);

mesh.r = [mesh.r Bore.rad]; %mesh lines in radial direction
mesh.z = [-Bore.length/2 mesh.z Bore.length/2]; %mesh lines in z-direction

mesh = SmoothMesh(mesh, mesh_res, 1.5);

%% check the cell limit
numCells = numel(mesh.r)*numel(mesh.a)*numel(mesh.z);

%% define human body model
CSX = AddDiscMaterial(CSX, 'body_model', 'File', body_model_file, 'Scale', 1/unit, 'Transform', body_model_transform);
start = [mesh.r(1)   mesh.a(1)   mesh.z(1)];
stop =  [mesh.r(end) mesh.a(end) mesh.z(end)];
CSX = AddBox(CSX, 'body_model', 0, start, stop);


%% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start = [0      mesh.a(1)   -BC.length/2];
stop =  [BC.rad mesh.a(end) +BC.length/2];

CSX = AddDump(CSX,'Ef','FileType',1,'DumpType',10,'DumpMode',2,'Frequency',f0);
CSX = AddBox(CSX,'Ef',0 , start,stop);

CSX = AddDump(CSX,'Hf','FileType',1,'DumpType',11,'DumpMode',2,'Frequency',f0);
CSX = AddBox(CSX,'Hf',0 , start,stop);

CSX = AddDump(CSX,'SAR','FileType',1,'DumpType',20,'DumpMode',2,'Frequency',f0);
CSX = AddBox(CSX,'SAR',0 , start,stop);

start = [0      mesh.a(1)   0];
stop =  [BC.rad mesh.a(end) 0];
CSX = AddDump(CSX,'Ht','FileType',1,'DumpType',1,'DumpMode',2);
CSX = AddBox(CSX,'Ht',0 , start,stop);

%% finalize mesh
% add some lines for the pml in +/- z- direction
mesh = AddPML(mesh, [0 0 0 0 10 10], 1);

% define the mesh
CSX = DefineRectGrid(CSX, unit, mesh);

%% Write file & run openEMS
Sim_Path = ['tmp_' mfilename];

if (postproc_only==0)
    [status, message, messageid] = rmdir(Sim_Path,'s'); %delete old results
    [status, message, messageid] = mkdir(Sim_Path);     %create folder

    WriteOpenEMS([Sim_Path '/BirdCage.xml'],FDTD,CSX);
end

if (GeomPlot==1)
    CSXGeomPlot( [Sim_Path '/BirdCage.xml'] , ['--export-polydata-vtk=' Sim_Path ' --RenderDiscMaterial -v']);
end

if (postproc_only==0)
    RunOpenEMS(Sim_Path, 'BirdCage.xml');
end

%%
freq = linspace(excite.f_0-excite.f_c,excite.f_0+excite.f_c,201);
port = calcPort(port, Sim_Path, freq);

close all
s11 = port{1}.uf.ref./port{1}.uf.inc;
s22 = port{2}.uf.ref./port{2}.uf.inc;

% the s-parameter may be larger than 1 (0dB) since all ports are excited
% and do not have a perfect port isolation
plot(freq*1e-6,20*log10(abs(s11)),'Linewidth',2)
hold on
grid on
plot(freq*1e-6,20*log10(abs(s22)),'r--','Linewidth',2)
legend('s11','s22');

%% read SAR values on a xy-plane (range)
[SAR SAR_mesh] = ReadHDF5Dump([Sim_Path '/SAR.h5'],'Range',{[],[],0},'CloseAlpha',1);
SAR = SAR.FD.values{1};

% SAR plot
figure()
[R A] = ndgrid(SAR_mesh.lines{1},SAR_mesh.lines{2});
X = R.*cos(A);Y = R.*sin(A);
colormap('hot');
h = pcolor(X,Y,(squeeze(SAR)));
% h = pcolor(X,Y,log10(squeeze(SAR)));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('local SAR');
axis equal tight

%% plot B1+/- on an xy-plane
[H_field H_mesh] = ReadHDF5Dump([Sim_Path '/Hf.h5'],'Range',{[0 0.1],[],0},'CloseAlpha',1);
% create a 2D grid to plot on
[R A] = ndgrid(H_mesh.lines{1},H_mesh.lines{2});
X = R.*cos(A);
Y = R.*sin(A);

% calc Bx,By (from Br and Ba), B1p, B1m
Bx = MUE0*(H_field.FD.values{1}(:,:,:,1).*cos(A) - H_field.FD.values{1}(:,:,:,2).*sin(A));
By = MUE0*(H_field.FD.values{1}(:,:,:,1).*sin(A) + H_field.FD.values{1}(:,:,:,2).*cos(A));
B1p = 0.5*(Bx+1j*By);
B1m = 0.5*(Bx-1j*By);

Dump2VTK([Sim_Path '/B1p_xy.vtk'], abs(B1p), H_mesh, 'B-Field');
Dump2VTK([Sim_Path '/B1m_xy.vtk'], abs(B1m), H_mesh, 'B-Field');

maxB1 = max([abs(B1p(:)); abs(B1m(:))]);

% B1+ plot
figure()
subplot(1,2,1);
h = pcolor(X,Y,abs(B1p));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^+ field (dB)');
caxis([0 maxB1]);
axis equal tight

% B1- plot
subplot(1,2,2);
h = pcolor(X,Y,abs(B1m));
set(h,'EdgeColor','none');
xlabel('x -->');
ylabel('y -->');
title('B_1^- field (dB)');
caxis([0 maxB1]);
axis equal tight

%%
ConvertHDF5_VTK([Sim_Path '/Hf.h5'],[Sim_Path '/Hf_xy'],'Range',{[],[],0},'CloseAlpha',1)
ConvertHDF5_VTK([Sim_Path '/SAR.h5'],[Sim_Path '/SAR_xy'],'Range',{[],[],0},'CloseAlpha',1)
