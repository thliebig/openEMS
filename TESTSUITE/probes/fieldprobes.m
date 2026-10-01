function pass = fieldprobes(varargin)
%pass = fieldprobes(<key>, <value>, ...)
%
% Infinitesimal dipole in free space. E- and H-field point probes are compared
% against the same position taken out of an HDF5 field dump: both sample the
% engine through Engine_Interface_Base, so they have to agree to within the
% file precision.
%
% The probed amplitude is checked as well -- identical zeros would otherwise
% compare just fine.
%
% See ts_options for the accepted options.
%
% openEMS testsuite
% -----------------
%
% See also run_testsuite, ts_options

% the shared check helpers, also when this test is run on its own
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'helperscripts'));

opt = ts_options(varargin{:});
pass = 1;

physical_constants;

% LIMITS
limit_max_time_diff  = 1e-13;
limit_max_amp_diff   = 1e-7;  % relative amplitude difference
limit_min_e_amp      = 5e-3;
limit_min_h_amp      = 1e-7;


% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drawingunit = 1e-6; % specify everything in um
Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX = 'fieldprobes.xml';

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

% dump boxes on the six faces of a box around the dipole
s1 = [-4.5, -4.5, -4.5] * dipole_length/2;
s2 = [ 4.5,  4.5,  4.5] * dipole_length/2;
face_start = {s1, [s2(1) s1(2) s1(3)], s1, [s1(1) s2(2) s1(3)], s1, [s1(1) s1(2) s2(3)]};
face_stop  = {[s1(1) s2(2) s2(3)], s2, [s2(1) s1(2) s2(3)], s2, [s2(1) s2(2) s1(3)], s2};
face_name  = {'xn','xp','yn','yp','zn','zp'};

% the probe position on each face
coords = {[s1(1) 0 0], [s2(1) 0 0], [0 s1(2) 0], [0 s2(2) 0], [0 0 s1(3)], [0 0 s2(3)]};

for n=1:numel(face_name)
    CSX = AddDump( CSX, ['Et_' face_name{n}], 'DumpType', 0, 'DumpMode', 0, 'FileType', 1 );
    CSX = AddBox( CSX, ['Et_' face_name{n}], 0, face_start{n}, face_stop{n} );
    CSX = AddDump( CSX, ['Ht_' face_name{n}], 'DumpType', 1, 'DumpMode', 0, 'FileType', 1 );
    CSX = AddBox( CSX, ['Ht_' face_name{n}], 0, face_start{n}, face_stop{n} );

    CSX = AddProbe( CSX, ['et' num2str(n)], 2 );
    CSX = AddPoint( CSX, ['et' num2str(n)], 0, coords{n} );
    CSX = AddProbe( CSX, ['ht' num2str(n)], 3 );
    CSX = AddPoint( CSX, ['ht' num2str(n)], 0, coords{n} );
end


% setup FDTD parameters & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%
FDTD = InitFDTD( 'NrTS', 10000, 'EndCriteria', 1e-6, 'OverSampling', 10 );
FDTD = SetGaussExcite( FDTD, 0, f_max );
FDTD = SetBoundaryCond( FDTD, {'MUR','MUR','MUR','MUR','MUR','MUR'} );

% Write openEMS compatible xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

% run openEMS
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent = opt.Silent;
RunOpenEMS( Sim_Path, Sim_CSX, opt.openEMS_opts, Settings );


%% POSTPROCESS
% Both field types are checked the same way: pick the probe position out of
% the dump and compare it against the point probe of that position.
checks = { ...
    struct('name','E', 'prefix','Et_', 'probe','et', 'min_amp',limit_min_e_amp, 'amp_comp',3), ...
    struct('name','H', 'prefix','Ht_', 'probe','ht', 'min_amp',limit_min_h_amp, 'amp_comp',[1 2]) };

for c=1:numel(checks)
    chk = checks{c};
    len_ok  = 1;
    time_ok = 1;
    amp_ok  = 1;
    min_ok  = 1;
    max_time_diff = 0;
    max_amp_diff  = 0;
    min_amp = Inf;

    for n=1:numel(face_name)
        dump  = ReadHDF5FieldData( [Sim_Path '/' chk.prefix face_name{n} '.h5'] );
        dmesh = ReadHDF5Mesh(      [Sim_Path '/' chk.prefix face_name{n} '.h5'] );
        probe = load( [Sim_Path '/' chk.probe num2str(n)] );

        idx = [1 1 1];
        for d=1:3
            if numel(dmesh.lines{d}) > 1
                idx(d) = interp1( dmesh.lines{d}, 1:numel(dmesh.lines{d}), coords{n}(d), 'nearest' );
            end
        end
        if opt.Verbose
            fprintf('    %s %s: dump position (%g,%g,%g) m, indices (%d,%d,%d)\n', ...
                    chk.name, face_name{n}, dmesh.lines{1}(idx(1)), ...
                    dmesh.lines{2}(idx(2)), dmesh.lines{3}(idx(3)), idx(1), idx(2), idx(3));
        end

        field = zeros(numel(dump.TD.values),3);
        for t=1:numel(dump.TD.values)
            field(t,:) = squeeze(dump.TD.values{t}(idx(1),idx(2),idx(3),:));
        end
        field_t   = reshape( dump.TD.time, [], 1 );
        probe_val = probe(:,2:4);

        % the dump and the probe must cover the same time steps
        if size(field,1) ~= size(probe,1)
            len_ok = 0;
            break
        end
        max_time_diff = max( max_time_diff, max(abs(field_t - probe(:,1))) );

        % relative deviation where the dump is non-zero; a zero dump sample
        % must be matched by a zero probe sample
        rel_diff = zeros(size(field));
        nz = (field ~= 0);
        rel_diff(nz) = (field(nz) - probe_val(nz)) ./ field(nz);
        if any( probe_val(~nz) ~= 0 )
            rel_diff(~nz) = Inf;
        end
        max_amp_diff = max( max_amp_diff, max(abs(rel_diff(:))) );

        % peak amplitude per component, the weakest of them has to show up
        min_amp = min( min_amp, min(max(abs(field(:,chk.amp_comp)), [], 1)) );

        if opt.Plots
            figure
            for d=1:3
                subplot(2,3,d);
                plot( field_t, [field(:,d) probe(:,1+d)] );
                title([chk.name '_' char('x'+d-1) ' ' face_name{n}]);
                subplot(2,3,3+d);
                plot( field_t, rel_diff(:,d) );
            end
        end
    end

    time_ok = max_time_diff <= limit_max_time_diff;
    amp_ok  = max_amp_diff  <= limit_max_amp_diff;
    min_ok  = min_amp       >= chk.min_amp;

    pass = ts_check([chk.name '-probe/dump sample count'], len_ok) && pass;
    if ~len_ok
        continue
    end
    pass = ts_check([chk.name '-probe/dump time base'], time_ok, ...
                    'max. difference %.3g s (limit %.3g s)', max_time_diff, limit_max_time_diff) && pass;
    pass = ts_check([chk.name '-probe/dump amplitude'], amp_ok, ...
                    'max. rel. deviation %.3g (limit %.3g)', max_amp_diff, limit_max_amp_diff) && pass;
    pass = ts_check([chk.name '-field actually excited'], min_ok, ...
                    'weakest face amplitude %.3g (limit %.3g)', min_amp, chk.min_amp) && pass;
end

pass = ts_finish(opt, mfilename, pass, Sim_Path);
