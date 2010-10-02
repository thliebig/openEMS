function pass = cavity
%pass = cavity
% 
% Checks, if different engines produces identical results

CLEANUP = 1;        % if enabled and result is PASS, remove simulation folder
STOP_IF_FAILED = 1; % if enabled and result is FAILED, stop with error

engines = {'' '--engine=sse' '--engine=sse-compressed' '--engine=multithreaded'};

isOctave = exist('OCTAVE_VERSION','builtin') ~= 0;
if isOctave
    old_crr = confirm_recursive_rmdir(0);
end

global Sim_Path Sim_CSX
Sim_Path = 'tmp';
Sim_CSX = 'cavity.xml';

for n=1:numel(engines)
    result{n} = sim( engines{n} );
end

pass = compare( result );

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

if isOctave
    confirm_recursive_rmdir(old_crr);
end
return


function result = sim( openEMS_options )
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
FDTD = InitFDTD( 2000,1e-6 );
FDTD = SetGaussExcite(FDTD,(f_stop-f_start)/2,(f_stop-f_start)/2);
BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PEC' 'PEC'}; % PEC boundaries
FDTD = SetBoundaryCond(FDTD,BC);

% setup CSXCAD geometry
CSX = InitCSX();
mesh.x = linspace(0,a,26);
mesh.y = linspace(0,b,11);
mesh.z = linspace(0,d,32);
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
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_options );

% collect result
E.mesh = ReadHDF5Mesh( [Sim_Path '/Et.h5'] );
E.data = ReadHDF5FieldData( [Sim_Path '/Et.h5'] );
H.mesh = ReadHDF5Mesh( [Sim_Path '/Ht.h5'] );
H.data = ReadHDF5FieldData( [Sim_Path '/Ht.h5'] );
result.E = E;
result.H = H;




function pass = compare( results )
pass = 0;
% n=1: reference simulation
for n=2:numel(results)
    % iterate over all simulations
    EHfields = fieldnames(results{1});
    for m=1:numel(EHfields)
        % iterate over all fields (E, H)
        EHfield = EHfields{m};
        for o=1:numel(results{1}.(EHfield).data.values)
            % iterate over all timesteps
            if results{1}.(EHfield).data.values{o} ~= results{n}.(EHfield).data.values{o}
                disp( ['compare error: n=' num2str(n) '  field=' EHfield '  timestep:' num2str(o) '=' results{1}.(EHfield).data.names{o}] );
                disp( '  coords:' );
                find( results{1}.(EHfield).data.values{o} ~= results{n}.(EHfield).data.values{o} )
                return
            end
        end
    end
    disp( ['simulation ' num2str(n) ' is identical to simulation 1'] );
end
pass = 1;
