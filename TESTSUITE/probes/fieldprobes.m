function pass = fieldprobes( openEMS_options, options )
%
% infinitesimal dipole in free-space
%
% E/H-field probes are compared to hdf5 field dumps
%

pass = 1;

physical_constants;


ENABLE_PLOTS = 1;
CLEANUP = 1;        % if enabled and result is PASS, remove simulation folder
STOP_IF_FAILED = 1; % if enabled and result is FAILED, stop with error
VERBOSE = 1;
SILENT = 0;         % 0=show openEMS output

if nargin < 1
    openEMS_options = '';
end
if nargin < 2
    options = '';
end
if any(strcmp( options, 'run_testsuite' ))
    ENABLE_PLOTS = 0;
    STOP_IF_FAILED = 0;
    SILENT = 1;
    VERBOSE = 0;
end

% LIMITS
limit_max_time_diff  = 1e-13;
limit_max_amp_diff   = 1e-7;  %relative amplitude difference
limit_min_e_amp      = 5e-3;
limit_min_h_amp      = 1e-7;


% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawingunit = 1e-6; % specify everything in um
Sim_Path = 'tmp_fieldprobes';
Sim_CSX = 'tmp.xml';

f_max = 1e9;
lambda = c0/f_max /drawingunit;

% setup geometry values
dipole_length = lambda/50;


% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSX = InitCSX();
mesh.x = -dipole_length*20:dipole_length/2:dipole_length*20;
mesh.y = -dipole_length*20:dipole_length/2:dipole_length*20;
mesh.z = -dipole_length*20:dipole_length/2:dipole_length*20;
CSX = DefineRectGrid( CSX, drawingunit, mesh );

% excitation
CSX = AddExcitation( CSX, 'infDipole', 1, [0 0 1] );
start = [0, 0, -dipole_length/2];
stop  = [0, 0, +dipole_length/2];
CSX = AddBox( CSX, 'infDipole', 1, start, stop );

% NFFF contour
s1 = [-4.5, -4.5, -4.5] * dipole_length/2;
s2 = [ 4.5,  4.5,  4.5] * dipole_length/2;
CSX = AddBox( AddDump(CSX,'Et_xn','DumpType',0,'DumpMode',0,'FileType',1), 'Et_xn', 0, s1, [s1(1) s2(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Et_xp','DumpType',0,'DumpMode',0,'FileType',1), 'Et_xp', 0, [s2(1) s1(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Et_yn','DumpType',0,'DumpMode',0,'FileType',1), 'Et_yn', 0, s1, [s2(1) s1(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Et_yp','DumpType',0,'DumpMode',0,'FileType',1), 'Et_yp', 0, [s1(1) s2(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Et_zn','DumpType',0,'DumpMode',0,'FileType',1), 'Et_zn', 0, s1, [s2(1) s2(2) s1(3)] );
CSX = AddBox( AddDump(CSX,'Et_zp','DumpType',0,'DumpMode',0,'FileType',1), 'Et_zp', 0, [s1(1) s1(2) s2(3)], s2 );
CSX = AddBox( AddDump(CSX,'Ht_xn','DumpType',1,'DumpMode',0,'FileType',1), 'Ht_xn', 0, s1, [s1(1) s2(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Ht_xp','DumpType',1,'DumpMode',0,'FileType',1), 'Ht_xp', 0, [s2(1) s1(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Ht_yn','DumpType',1,'DumpMode',0,'FileType',1), 'Ht_yn', 0, s1, [s2(1) s1(2) s2(3)] );
CSX = AddBox( AddDump(CSX,'Ht_yp','DumpType',1,'DumpMode',0,'FileType',1), 'Ht_yp', 0, [s1(1) s2(2) s1(3)], s2 );
CSX = AddBox( AddDump(CSX,'Ht_zn','DumpType',1,'DumpMode',0,'FileType',1), 'Ht_zn', 0, s1, [s2(1) s2(2) s1(3)] );
CSX = AddBox( AddDump(CSX,'Ht_zp','DumpType',1,'DumpMode',0,'FileType',1), 'Ht_zp', 0, [s1(1) s1(2) s2(3)], s2 );

% E-field probes
coords{1} = [s1(1) 0 0];
CSX = AddPoint( AddProbe(CSX,'et1',2), 'et1', 0, coords{1} );
coords{2} = [s2(1) 0 0];
CSX = AddPoint( AddProbe(CSX,'et2',2), 'et2', 0, coords{2} );
coords{3} = [0 s1(2) 0];
CSX = AddPoint( AddProbe(CSX,'et3',2), 'et3', 0, coords{3} );
coords{4} = [0 s2(2) 0];
CSX = AddPoint( AddProbe(CSX,'et4',2), 'et4', 0, coords{4} );
coords{5} = [0 0 s1(3)];
CSX = AddPoint( AddProbe(CSX,'et5',2), 'et5', 0, coords{5} );
coords{6} = [0 0 s2(3)];
CSX = AddPoint( AddProbe(CSX,'et6',2), 'et6', 0, coords{6} );

% H-field probes
CSX = AddPoint( AddProbe(CSX,'ht1',3), 'ht1', 0, [s1(1) 0 0] );
CSX = AddPoint( AddProbe(CSX,'ht2',3), 'ht2', 0, [s2(1) 0 0] );
CSX = AddPoint( AddProbe(CSX,'ht3',3), 'ht3', 0, [0 s1(2) 0] );
CSX = AddPoint( AddProbe(CSX,'ht4',3), 'ht4', 0, [0 s2(2) 0] );
CSX = AddPoint( AddProbe(CSX,'ht5',3), 'ht5', 0, [0 0 s1(3)] );
CSX = AddPoint( AddProbe(CSX,'ht6',3), 'ht6', 0, [0 0 s2(3)] );



% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_timesteps = 10000;
min_decrement = 1e-6;
FDTD = InitFDTD( max_timesteps, min_decrement,'OverSampling',10 );
FDTD = SetGaussExcite( FDTD, 0, f_max );
BC = [2 2 2 2 2 2];
FDTD = SetBoundaryCond( FDTD, BC );

% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~] = rmdir(Sim_Path,'s');
[~,~,~] = mkdir(Sim_Path);
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

% run openEMS
folder = fileparts( mfilename('fullpath') );
Settings.LogFile = [folder '/' Sim_Path '/openEMS.log'];
Settings.Silent = SILENT;
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_options, Settings );


%% POSTPROCESS
filenames_E = {'Et_xn.h5','Et_xp.h5','Et_yn.h5','Et_yp.h5','Et_zn.h5','Et_zp.h5'};
filenames_H = {'Ht_xn.h5','Ht_xp.h5','Ht_yn.h5','Ht_yp.h5','Ht_zn.h5','Ht_zp.h5'};

for n=1:numel(filenames_E)
    Et{n} = ReadHDF5FieldData( [Sim_Path '/' filenames_E{n}] );
    E_mesh{n} = ReadHDF5Mesh( [Sim_Path '/' filenames_E{n}] );
    Ht{n} = ReadHDF5FieldData( [Sim_Path '/' filenames_H{n}] );
    H_mesh{n} = ReadHDF5Mesh( [Sim_Path '/' filenames_H{n}] );
    Et_probe{n} = load( [Sim_Path '/et' num2str(n)] );
    Ht_probe{n} = load( [Sim_Path '/ht' num2str(n)] );
end

if ENABLE_PLOTS
    close all
end

%
% E-fields
%
if VERBOSE, disp( 'extracting field components from field dumps...' ); end
for n=1:6
    if numel(E_mesh{n}.lines{1}) > 1
        x_idx = interp1( E_mesh{n}.lines{1}, 1:numel(E_mesh{n}.lines{1}), coords{n}(1), 'nearest' );
    else
        x_idx = 1;
    end
    if numel(E_mesh{n}.lines{2}) > 1
        y_idx = interp1( E_mesh{n}.lines{2}, 1:numel(E_mesh{n}.lines{2}), coords{n}(2), 'nearest' );
    else
        y_idx = 1;
    end
    if numel(E_mesh{n}.lines{3}) > 1
        z_idx = interp1( E_mesh{n}.lines{3}, 1:numel(E_mesh{n}.lines{3}), coords{n}(3), 'nearest' );
    else
        z_idx = 1;
    end
    
    if VERBOSE
        disp( ['n=' num2str(n) ' coords: (' num2str(E_mesh{n}.lines{1}(x_idx)) ','...
               num2str(E_mesh{n}.lines{2}(y_idx)) ','...
               num2str(E_mesh{n}.lines{3}(z_idx)) ') m  indices: ('...
               num2str(x_idx) ',' num2str(y_idx) ',' num2str(z_idx) ')'] );
    end
    
    field_x = zeros(numel(Et{n}.TD.values),1);
    field_y = zeros(numel(Et{n}.TD.values),1);
    field_z = zeros(numel(Et{n}.TD.values),1);
    for t=1:numel(Et{n}.TD.values)
        field_x(t) = squeeze(Et{n}.TD.values{t}(x_idx,y_idx,z_idx,1));
        field_y(t) = squeeze(Et{n}.TD.values{t}(x_idx,y_idx,z_idx,2));
        field_z(t) = squeeze(Et{n}.TD.values{t}(x_idx,y_idx,z_idx,3));
    end
    field_t = reshape( Et{n}.TD.time, [], 1 );
    
    % check vector length
    if numel(field_x) ~= size(Et_probe{n},1)
        pass = 0;
        disp( 'probes/fieldprobes.m (vector length):  * FAILED *' );
        break
    end
    
    % check absolute simulation time
    if any(abs(field_t - Et_probe{n}(:,1)) > limit_max_time_diff)
        pass = 0;
        disp( 'probes/fieldprobes.m (time inconsistant):  * FAILED *' );
        break
    end
    
    if ENABLE_PLOTS
        figure
        subplot(2,3,1);
        plot( field_t, [field_x Et_probe{n}(:,2)] );
        subplot(2,3,2);
        plot( field_t, [field_y Et_probe{n}(:,3)] );
        subplot(2,3,3);
        plot( field_t, [field_z Et_probe{n}(:,4)] );
        subplot(2,3,4);
        plot( field_t, (field_x - Et_probe{n}(:,2))./field_x );
        subplot(2,3,5);
        plot( field_t, (field_y - Et_probe{n}(:,3))./field_y );
        subplot(2,3,6);
        plot( field_t, (field_z - Et_probe{n}(:,4))./field_z );
    end
    
    % difference
    if any( abs( (field_x - Et_probe{n}(:,2))./field_x) > limit_max_amp_diff ) || ...
       any( abs( (field_y - Et_probe{n}(:,3))./field_y) > limit_max_amp_diff ) || ...
       any( abs( (field_z - Et_probe{n}(:,4))./field_z) > limit_max_amp_diff )
        pass = 0;
        disp( 'probes/fieldprobes.m (amplitudes differ too much):  * FAILED *' );
        break
    end
    
    % check absolute field strength of z component
    if max(abs(field_z)) < limit_min_e_amp
        pass = 0;
        disp( 'probes/fieldprobes.m (amplitude of z-component too small):  * FAILED *' );
        break
    end
end

%
% H-fields
%
if VERBOSE, disp( 'extracting field components from field dumps...' ); end
for n=1:6
    if numel(H_mesh{n}.lines{1}) > 1
        x_idx = interp1( H_mesh{n}.lines{1}, 1:numel(H_mesh{n}.lines{1}), coords{n}(1), 'nearest' );
    else
        x_idx = 1;
    end
    if numel(E_mesh{n}.lines{2}) > 1
        y_idx = interp1( H_mesh{n}.lines{2}, 1:numel(H_mesh{n}.lines{2}), coords{n}(2), 'nearest' );
    else
        y_idx = 1;
    end
    if numel(E_mesh{n}.lines{3}) > 1
        z_idx = interp1( H_mesh{n}.lines{3}, 1:numel(H_mesh{n}.lines{3}), coords{n}(3), 'nearest' );
    else
        z_idx = 1;
    end
    
    if VERBOSE
        disp( ['n=' num2str(n) ' coords: (' num2str(E_mesh{n}.lines{1}(x_idx)) ','...
               num2str(E_mesh{n}.lines{2}(y_idx)) ','...
               num2str(E_mesh{n}.lines{3}(z_idx)) ') m  indices: ('...
               num2str(x_idx) ',' num2str(y_idx) ',' num2str(z_idx) ')'] );
    end
    
    field_x = zeros(numel(Ht{n}.TD.values),1);
    field_y = zeros(numel(Ht{n}.TD.values),1);
    field_z = zeros(numel(Ht{n}.TD.values),1);
    for t=1:numel(Ht{n}.TD.values)
        field_x(t) = squeeze(Ht{n}.TD.values{t}(x_idx,y_idx,z_idx,1));
        field_y(t) = squeeze(Ht{n}.TD.values{t}(x_idx,y_idx,z_idx,2));
        field_z(t) = squeeze(Ht{n}.TD.values{t}(x_idx,y_idx,z_idx,3));
    end
    field_t = reshape( Ht{n}.TD.time, [], 1 );
    
    % check vector length
    if numel(field_x) ~= size(Ht_probe{n},1)
        pass = 0;
        disp( 'probes/fieldprobes.m (vector length):  * FAILED *' );
        break
    end
    
    % check absolute simulation time
    if any(abs(field_t - Ht_probe{n}(:,1)) > limit_max_time_diff)
        pass = 0;
        disp( 'probes/fieldprobes.m (time inconsistant):  * FAILED *' );
        break
    end
    
    if ENABLE_PLOTS
        figure
        subplot(2,3,1);
        plot( field_t, [field_x Ht_probe{n}(:,2)] );
        subplot(2,3,2);
        plot( field_t, [field_y Ht_probe{n}(:,3)] );
        subplot(2,3,3);
        plot( field_t, [field_z Ht_probe{n}(:,4)] );
        subplot(2,3,4);
        plot( field_t, (field_x - Ht_probe{n}(:,2))./field_x );
        subplot(2,3,5);
        plot( field_t, (field_y - Ht_probe{n}(:,3))./field_y );
        subplot(2,3,6);
        plot( field_t, (field_z - Ht_probe{n}(:,4))./field_z );
    end
    
    % difference
    if any( abs( (field_x - Ht_probe{n}(:,2))./field_x) > limit_max_amp_diff ) || ...
       any( abs( (field_y - Ht_probe{n}(:,3))./field_y) > limit_max_amp_diff ) || ...
       any( abs( (field_z - Ht_probe{n}(:,4))./field_z) > limit_max_amp_diff )
        pass = 0;
        disp( 'probes/fieldprobes.m (amplitudes differ too much):  * FAILED *' );
        break
    end
    
    % check absolute field strength of z component
    if (max(abs(field_x)) < limit_min_h_amp) || (max(abs(field_y)) < limit_min_h_amp)
        pass = 0;
        disp( 'probes/fieldprobes.m (amplitude of x- or y-component too small):  * FAILED *' );
        break
    end
end




if pass
    disp( 'probes/fieldprobes.m:  pass' );
end

    
if pass && CLEANUP
    rmdir( Sim_Path, 's' );
end
if ~pass && STOP_IF_FAILED
    error 'test failed';
end
