function pass = cavity( openEMS_options, options )
%pass = cavity( openEMS_options, options )
% 
% Checks, if different engines produces identical results

CLEANUP = 1;        % if enabled and result is PASS, remove simulation folder
STOP_IF_FAILED = 1; % if enabled and result is FAILED, stop with error
global ENABLE_PLOTS;
ENABLE_PLOTS = 1;
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
end
% clean openEMS_options
openEMS_options = regexprep( openEMS_options, '--engine=\w+', '' );

engines = {'--engine=basic' '--engine=sse' '--engine=sse-compressed' '--engine=multithreaded'};
% engines = [engines {'--engine=sse-compressed-linear' '--engine=multithreaded-linear'}];

global Sim_Path Sim_CSX
Sim_Path = 'tmp_cavity';
Sim_CSX = 'cavity.xml';

for n=1:numel(engines)
    result{n} = sim( [engines{n} ' ' openEMS_options], SILENT );
end

pass = compare( result, SILENT );

if pass
    disp( 'enginetests/cavity.m (engine comparison):  pass' );
else
    disp( 'enginetests/cavity.m (engine comparison):  * FAILED *' );
end

if pass && CLEANUP
    rmdir( Sim_Path, 's' );
end
if ~pass && STOP_IF_FAILED
    error 'test failed'
end

return


function result = sim( openEMS_options, SILENT )
global Sim_Path Sim_CSX
physical_constants;

% structure
a = 5e-2;
b = 2e-2;
d = 6e-2;
if ~((b<a) && (a<d))
    error 'correct the dimensions of the cavity'
end

f_start = 1e9;
f_stop = 10e9;

% prepare simulation dir
[status,message,messageid] = rmdir(Sim_Path,'s');
[status,message,messageid] = mkdir(Sim_Path);

% setup FDTD parameter
FDTD = InitFDTD( 1000, 0 );
FDTD = SetGaussExcite(FDTD,(f_stop-f_start)/2,(f_stop-f_start)/2);
BC = {'MUR' 'PML_8' 'PMC' 'PEC' 'PEC' 'PEC'}; % boundaries
FDTD = SetBoundaryCond(FDTD,BC);

% setup CSXCAD geometry
CSX = InitCSX();
mesh.x = linspace(0,a,27);
mesh.y = linspace(0,b,11);
mesh.z = linspace(0,d,33);
CSX = DefineRectGrid(CSX, 1,mesh);

% excitation
CSX = AddExcitation(CSX,'excite1',0,[1 1 1]);
p(1,1) = mesh.x(floor(end*2/3));
p(2,1) = mesh.y(floor(end*2/3));
p(3,1) = mesh.z(floor(end*2/3));
p(1,2) = mesh.x(floor(end*2/3)+1);
p(2,2) = mesh.y(floor(end*2/3)+1);
p(3,2) = mesh.z(floor(end*2/3)+1);
CSX = AddCurve( CSX, 'excite1', 0, p );

% probes
CSX = AddProbe( CSX, 'E_probe', 2 );
p(1,1) = mesh.x(floor(end*1/3));
p(2,1) = mesh.y(floor(end*1/3));
p(3,1) = mesh.z(floor(end*1/3));
CSX = AddPoint( CSX, 'E_probe', 0, p );
CSX = AddProbe( CSX, 'H_probe', 3 );
CSX = AddPoint( CSX, 'H_probe', 0, p );

% material
CSX = AddMaterial( CSX, 'RO4350B', 'Epsilon', 3.66 );
start = [mesh.x(3) mesh.y(3) mesh.z(3)];
stop  = [mesh.x(5) mesh.y(4) mesh.z(6)];
CSX = AddBox( CSX, 'RO4350B', 100, start, stop );

% dump
CSX = AddDump( CSX, 'Et', 'DumpType', 0, 'DumpMode', 0, 'FileType', 1 ); % hdf5 E-field dump without interpolation
pos1 = [mesh.x(1) mesh.y(1) mesh.z(1)];
pos2 = [mesh.x(end) mesh.y(end) mesh.z(end)];
CSX = AddBox( CSX, 'Et', 0, pos1, pos2 );

% dump
CSX = AddDump( CSX, 'Ht', 'DumpType', 1, 'DumpMode', 0, 'FileType', 1 ); % hdf5 H-field dump without interpolation
pos1 = [mesh.x(1) mesh.y(1) mesh.z(1)];       % should be half a cell more than now
pos2 = [mesh.x(end) mesh.y(end) mesh.z(end)]; % should be half a cell less than now
CSX = AddBox( CSX, 'Ht', 0, pos1, pos2 );

% Write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% cd to working dir and run openEMS
folder = fileparts( mfilename('fullpath') );
Settings.LogFile = [folder '/' Sim_Path '/openEMS.log'];
Settings.Silent = SILENT;
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_options, Settings );

% collect result
E.mesh = ReadHDF5Mesh( [Sim_Path '/Et.h5'] );
E.data = ReadHDF5FieldData( [Sim_Path '/Et.h5'] );
H.mesh = ReadHDF5Mesh( [Sim_Path '/Ht.h5'] );
H.data = ReadHDF5FieldData( [Sim_Path '/Ht.h5'] );
result.E = E;
result.H = H;
result.probes = ReadUI( {'E_probe','H_probe'}, Sim_Path );



function pass = compare( results, SILENT )
pass = 0;
% n=1: reference simulation
for n=2:numel(results)
    % iterate over all simulations
    EHfields = {'E','H'};
    for m=1:numel(EHfields)
        % iterate over all fields (E, H)
        EHfield = EHfields{m};
        for o=1:numel(results{1}.(EHfield).data.TD.values)
            % iterate over all timesteps
            cmp_result = results{1}.(EHfield).data.TD.values{o} ~= results{n}.(EHfield).data.TD.values{o};
            if any(cmp_result(:))
                disp( ['compare error: n=' num2str(n) '  field=' EHfield '  timestep:' num2str(o) '=' results{1}.(EHfield).data.names{o}] );
                disp( '  coords:' );
                find( results{1}.(EHfield).data.TD.values{o} ~= results{n}.(EHfield).data.TD.values{o} )
                return
            end
        end
    end
    if ~SILENT
        disp( ['simulation ' num2str(n) ' is identical to simulation 1'] );
    end
end

global ENABLE_PLOTS;
if ENABLE_PLOTS
    figure
    l = {};
    for n=1:numel(results)
        plot( results{n}.probes.TD{1}.t, results{n}.probes.TD{1}.val );
        hold all
        l = [l ['sim ' num2str(n)]];
    end
    legend( l );

    figure
    l = {};
    for n=1:numel(results)
        plot( results{n}.probes.TD{2}.t, results{n}.probes.TD{2}.val );
        hold all
        l = [l ['sim ' num2str(n)]];
    end
    legend( l );
end

pass = 1;
