function pass = cavity(varargin)
%pass = cavity(<key>, <value>, ...)
%
% Rectangular PEC cavity, excited by a short current curve. The spectrum of
% three orthogonal voltage probes is checked against the analytic resonances
%
%     f_mnl = c0/(2*pi) * sqrt((m*pi/a)^2 + (n*pi/b)^2 + (l*pi/d)^2)
%
% Two things are verified: every analytic resonance carries a peak ('inside'),
% and the spectrum stays low between them ('outside').
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
physical_constants;

% LIMITS - inside
lower_rel_limit = 1.3e-3;    % -0.13%
upper_rel_limit = 1.3e-3;    % +0.13%
lower_rel_limit_TM = 2.5e-3; % -0.25%
upper_rel_limit_TM = 0;      % +0%
min_rel_amplitude = 0.6;     % 60%
min_rel_amplitude_TM = 0.27; % 27%

% LIMITS - outside
outer_rel_limit = 0.02;
max_rel_amplitude = 0.17;


% structure
a = 5e-2;
b = 2e-2;
d = 6e-2;
if ~((b<a) && (a<d))
    error 'correct the dimensions of the cavity'
end

f_start = 1e9;
f_stop = 10e9;

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX = 'cavity.xml';

%setup FDTD parameter
FDTD = InitFDTD('NrTS', 20000, 'EndCriteria', 1e-6);
FDTD = SetGaussExcite(FDTD,(f_stop-f_start)/2,(f_stop-f_start)/2);
FDTD = SetBoundaryCond(FDTD,{'PEC','PEC','PEC','PEC','PEC','PEC'});

%setup CSXCAD geometry
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

%voltage calc
CSX = AddProbe(CSX,'ut1x',0);
pos1 = [mesh.x(floor(end/4)) mesh.y(floor(end/2))   mesh.z(floor(end/5))];
pos2 = [mesh.x(floor(end/4)+1) mesh.y(floor(end/2)) mesh.z(floor(end/5))];
CSX = AddBox(CSX,'ut1x', 0 ,pos1,pos2);

CSX = AddProbe(CSX,'ut1y',0);
pos1 = [mesh.x(floor(end/4)) mesh.y(floor(end/2))   mesh.z(floor(end/5))];
pos2 = [mesh.x(floor(end/4)) mesh.y(floor(end/2)+1) mesh.z(floor(end/5))];
CSX = AddBox(CSX,'ut1y', 0 ,pos1,pos2);

CSX = AddProbe(CSX,'ut1z',0);
pos1 = [mesh.x(floor(end/2)) mesh.y(floor(end/2)) mesh.z(floor(end/5))];
pos2 = [mesh.x(floor(end/2)) mesh.y(floor(end/2)) mesh.z(floor(end/5)+1)];
CSX = AddBox(CSX,'ut1z', 0 ,pos1,pos2);

%Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

% run openEMS
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent = opt.Silent;
RunOpenEMS( Sim_Path, Sim_CSX, opt.openEMS_opts, Settings );
UI = ReadUI( {[Sim_Path '/ut1x'], [Sim_Path '/ut1y'], [Sim_Path '/ut1z']} );



%
% analysis
%

% remove excitation from time series
t_start = 7e-10; % FIXME to be calculated
t_idx_start = interp1( UI.TD{1}.t, 1:numel(UI.TD{1}.t), t_start, 'nearest' );
for n=1:numel(UI.TD)
    UI.TD{n}.t = UI.TD{n}.t(t_idx_start:end);
    UI.TD{n}.val = UI.TD{n}.val(t_idx_start:end);
    [UI.FD{n}.f,UI.FD{n}.val] = FFT_time2freq( UI.TD{n}.t, UI.TD{n}.val );
end


f = UI.FD{1}.f;
ux = UI.FD{1}.val;
uy = UI.FD{2}.val;
uz = UI.FD{3}.val;

f_idx_start = interp1( f, 1:numel(f), f_start, 'nearest' );
f_idx_stop  = interp1( f, 1:numel(f), f_stop,  'nearest' );
f = f(f_idx_start:f_idx_stop);
ux = ux(f_idx_start:f_idx_stop);
uy = uy(f_idx_start:f_idx_stop);
uz = uz(f_idx_start:f_idx_stop);

% analytic formula for resonant wavenumber
k = @(m,n,l) sqrt( (m*pi/a)^2 + (n*pi/b)^2 + (l*pi/d)^2 );
f_TE101 = c0/(2*pi) * k(1,0,1);
f_TE102 = c0/(2*pi) * k(1,0,2);
f_TE201 = c0/(2*pi) * k(2,0,1);
f_TE202 = c0/(2*pi) * k(2,0,2);
f_TM110 = c0/(2*pi) * k(1,1,0);
f_TM111 = c0/(2*pi) * k(1,1,1);

f_TE = [f_TE101 f_TE102 f_TE201 f_TE202];
f_TM = [f_TM110 f_TM111];

% calculate frequency limits
temp = [f_start f_TE f_stop];
f_outer1 = [];
f_outer2 = [];
for n=1:numel(temp)-1
    f_outer1 = [f_outer1 temp(n) .* (1+outer_rel_limit)];
    f_outer2 = [f_outer2 temp(n+1) .* (1-outer_rel_limit)];
end

temp = [f_start f_TM f_stop];
f_outer1_TM = [];
f_outer2_TM = [];
for n=1:numel(temp)-1
    f_outer1_TM = [f_outer1_TM temp(n) .* (1+outer_rel_limit)];
    f_outer2_TM = [f_outer2_TM temp(n+1) .* (1-outer_rel_limit)];
end


if opt.Plots
    figure
    plot(f/1e9,abs(uy))
    max1 = max(abs(uy));
    hold on
    plot( repmat(f_TE,2,1)/1e9, repmat([0; max1],1,numel(f_TE)), 'm-.', 'LineWidth', 2 )
    plot( (repmat(f_TE,2,1) .* repmat(1-lower_rel_limit,2,numel(f_TE)))/1e9, repmat([0; max1],1,numel(f_TE)), 'r-', 'LineWidth', 1 )
    plot( (repmat(f_TE,2,1) .* repmat(1+upper_rel_limit,2,numel(f_TE)))/1e9, repmat([0; max1],1,numel(f_TE)), 'r-', 'LineWidth', 1 )
    plot( (repmat(f_TE,2,1) .* repmat([1-outer_rel_limit;1+outer_rel_limit],1,numel(f_TE)))/1e9, repmat(max1*min_rel_amplitude,2,numel(f_TE)), 'r-', 'LineWidth', 1 ) % freq limits
    plot( [f_outer1;f_outer2]/1e9, repmat(max1*max_rel_amplitude,2,numel(f_outer1)), 'g-', 'LineWidth', 1 ) % amplitude limits
    xlabel('Frequency (GHz)')
    legend( {'u_y','theoretical'} )
    title( 'TE-modes' )

    figure
    plot(f/1e9,abs(uz))
    max1 = max(abs(uz));
    hold on
    plot( repmat(f_TM,2,1)/1e9, repmat([0; max1],1,numel(f_TM)), 'm-.', 'LineWidth', 2 )
    plot( (repmat(f_TM,2,1) .* repmat(1-lower_rel_limit_TM,2,numel(f_TM)))/1e9, repmat([0; max1],1,numel(f_TM)), 'r-', 'LineWidth', 1 )
    plot( (repmat(f_TM,2,1) .* repmat(1+upper_rel_limit_TM,2,numel(f_TM)))/1e9, repmat([0; max1],1,numel(f_TM)), 'r-', 'LineWidth', 1 )
    plot( (repmat(f_TM,2,1) .* repmat([1-lower_rel_limit_TM;1+upper_rel_limit_TM],1,numel(f_TM)))/1e9, repmat(max1*min_rel_amplitude_TM,2,numel(f_TM)), 'r-', 'LineWidth', 1 ) % freq limits
    plot( [f_outer1_TM;f_outer2_TM]/1e9, repmat(max1*max_rel_amplitude,2,numel(f_outer1_TM)), 'g-', 'LineWidth', 1 ) % amplitude limits
    xlabel('Frequency (GHz)')
    legend( {'u_z','theoretical'} )
    title( 'TM-modes' )
end

pass = ts_check('TE resonances present', ...
         check_frequency( f, abs(uy), f_TE*(1+upper_rel_limit), f_TE*(1-lower_rel_limit), min_rel_amplitude, 'inside' ), ...
         '%s GHz, at least %g%% of the peak amplitude', mat2str(round(f_TE/1e8)/10), min_rel_amplitude*100);
pass = ts_check('TM resonances present', ...
         check_frequency( f, abs(uz), f_TM*(1+upper_rel_limit_TM), f_TM*(1-lower_rel_limit_TM), min_rel_amplitude_TM, 'inside' ), ...
         '%s GHz, at least %g%% of the peak amplitude', mat2str(round(f_TM/1e8)/10), min_rel_amplitude_TM*100) && pass;
pass = ts_check('no spurious TE response', ...
         check_frequency( f, abs(uy), f_outer2, f_outer1, max_rel_amplitude, 'outside' ), ...
         'below %g%% of the peak amplitude between the resonances', max_rel_amplitude*100) && pass;
pass = ts_check('no spurious TM response', ...
         check_frequency( f, abs(uz), f_outer2_TM, f_outer1_TM, max_rel_amplitude, 'outside' ), ...
         'below %g%% of the peak amplitude between the resonances', max_rel_amplitude*100) && pass;

pass = ts_finish(opt, mfilename, pass, Sim_Path);
