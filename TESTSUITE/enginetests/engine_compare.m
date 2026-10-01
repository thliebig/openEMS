function pass = engine_compare(varargin)
%pass = engine_compare(<key>, <value>, ...)
%
% Check that all engines produce bit-identical results.
%
% The same small cavity -- a dielectric block, mixed boundary conditions, E-
% and H-field dumps and point probes -- is run with every engine and the field
% dumps are compared element by element. The engines differ only in how they
% store and traverse the operator, so any difference at all is a bug.
%
% This test sweeps the engines itself and ignores the 'openEMS_opts' option.
% It also needs no '--exact-endcriteria': with 'EndCriteria', 0 every run is
% exactly NrTS timesteps long, which is both deterministic already and what
% makes the field dumps comparable at all. See ts_options for the rest.
%
% openEMS testsuite
% -----------------
%
% See also run_testsuite, ts_options

% the shared check helpers, also when this test is run on its own
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helperscripts'));

opt = ts_options(varargin{:});

engines = {'--engine=basic' '--engine=sse' '--engine=sse-compressed' ...
           '--engine=multithreaded'};

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX  = 'cavity.xml';

for n=1:numel(engines)
    if opt.Verbose
        fprintf('    running %s\n', engines{n});
    end
    result{n} = sim( Sim_Path, Sim_CSX, engines{n}, opt );
end

pass = compare( result, engines, opt );

pass = ts_finish(opt, mfilename, pass, Sim_Path);

end


function result = sim( Sim_Path, Sim_CSX, openEMS_options, opt )
physical_constants;

% structure
a = 5e-2;
b = 2e-2;
d = 6e-2;

f_start = 1e9;
f_stop = 10e9;

% setup FDTD parameter
FDTD = InitFDTD( 'NrTS', 1000, 'EndCriteria', 0 );
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
CleanupSimPath( Sim_Path );
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

% run openEMS
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent = opt.Silent;
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_options, Settings );

% collect result
E.mesh = ReadHDF5Mesh( [Sim_Path '/Et.h5'] );
E.data = ReadHDF5FieldData( [Sim_Path '/Et.h5'] );
H.mesh = ReadHDF5Mesh( [Sim_Path '/Ht.h5'] );
H.data = ReadHDF5FieldData( [Sim_Path '/Ht.h5'] );
result.E = E;
result.H = H;
result.probes = ReadUI( {'E_probe','H_probe'}, Sim_Path );

end


function pass = compare( results, engines, opt )
pass = 1;
% n=1: reference simulation
for n=2:numel(results)
    % iterate over all simulations
    EHfields = {'E','H'};
    identical = 1;
    detail = '';
    for m=1:numel(EHfields)
        % iterate over all fields (E, H)
        EHfield = EHfields{m};
        for o=1:numel(results{1}.(EHfield).data.TD.values)
            % iterate over all timesteps
            ref = results{1}.(EHfield).data.TD.values{o};
            cmp = results{n}.(EHfield).data.TD.values{o};
            if any(ref(:) ~= cmp(:))
                identical = 0;
                detail = sprintf('%s-field differs at timestep %s in %d of %d samples, max. difference %g', ...
                                 EHfield, results{1}.(EHfield).data.names{o}, ...
                                 sum(ref(:) ~= cmp(:)), numel(ref), ...
                                 max(abs(ref(:) - cmp(:))));
                break
            end
        end
        if ~identical, break; end
    end
    if identical
        pass = ts_check(sprintf('%s == %s', engines{n}, engines{1}), 1, ...
                        'all field dump samples identical') && pass;
    else
        pass = ts_check(sprintf('%s == %s', engines{n}, engines{1}), 0, ...
                        '%s', detail) && pass;
    end
end

if opt.Plots
    for p=1:2
        figure
        l = {};
        for n=1:numel(results)
            plot( results{n}.probes.TD{p}.t, results{n}.probes.TD{p}.val );
            hold all
            l = [l engines{n}];
        end
        legend( l );
    end
end

end
