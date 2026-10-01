function pass = nf2ff_dipole(varargin)
%pass = nf2ff_dipole(<key>, <value>, ...)
%
% Infinitesimal (Hertzian) z-dipole in free space, transformed to the far field
% through a near-field to far-field box. The far field of a Hertzian dipole is
% known in closed form, which makes this a check of the whole chain: the field
% dumps on the six faces, the nf2ff tool, and CreateNF2FFBox/CalcNF2FF.
%
%     E_theta ~ sin(theta),   E_phi = 0,   D_max = 1.5  (1.76 dBi)
%
% A face of the box picked up with the wrong sign or orientation breaks the
% pattern symmetry long before it breaks the directivity, so both are checked.
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

% LIMITS
limit_Dmax     = 0.03;   % 3% on the directivity of 1.5
limit_pattern  = 0.03;   % 3% on the normalized sin(theta) pattern
limit_symmetry = 0.02;   % 2% between two phi cuts of a rotationally
                         % symmetric pattern
limit_crosspol = 0.02;   % |E_phi| / max|E_theta|

% setup, in um
drawingunit   = 1e-6;
f_max         = 1e9;
lambda        = c0/f_max/drawingunit;
dipole_length = lambda/50;

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX  = 'nf2ff_dipole.xml';

%% geometry and mesh
CSX = InitCSX();
mesh.x = -dipole_length*10:dipole_length/2:dipole_length*10;
mesh.y = mesh.x;
mesh.z = mesh.x;

% the dipole: a hard z-directed current over one mesh line
CSX = AddExcitation( CSX, 'infDipole', 1, [0 0 1] );
start = [0 0 -dipole_length/2] - [0.1 0.1 0.1]*dipole_length/2;
stop  = [0 0 +dipole_length/2] + [0.1 0.1 0.1]*dipole_length/2;
CSX = AddBox( CSX, 'infDipole', 1, start, stop );

% the nf2ff box on the outermost mesh lines, PML behind it
[CSX, nf2ff] = CreateNF2FFBox( CSX, 'nf2ff', ...
                               [mesh.x(1)   mesh.y(1)   mesh.z(1)], ...
                               [mesh.x(end) mesh.y(end) mesh.z(end)] );
mesh = AddPML( mesh, 8 );
CSX = DefineRectGrid( CSX, drawingunit, mesh );

%% FDTD setup
FDTD = InitFDTD( 'NrTS', 2000, 'EndCriteria', 1e-6, 'OverSampling', 10 );
FDTD = SetGaussExcite( FDTD, f_max/2, f_max/2 );
FDTD = SetBoundaryCond( FDTD, repmat({'PML_8'}, 1, 6) );

%% run
WriteOpenEMS([Sim_Path '/' Sim_CSX], FDTD, CSX);
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent  = opt.Silent;
RunOpenEMS(Sim_Path, Sim_CSX, opt.openEMS_opts, Settings);

%% far field
theta = (0:2:180)/180*pi;
phi   = [0 pi/2];
nf2ff = CalcNF2FF( nf2ff, Sim_Path, f_max, theta, phi, 'Mode', 1, ...
                   'Verbose', 0 );

E_theta = nf2ff.E_theta{1};   % (theta, phi)
E_phi   = nf2ff.E_phi{1};

if opt.Plots
    figure
    plot(theta/pi*180, abs(E_theta)/max(abs(E_theta(:))), 'LineWidth', 2);
    hold on
    plot(theta/pi*180, sin(theta), 'k--', 'LineWidth', 2);
    xlabel('theta (deg)'); ylabel('|E_\theta| (normalized)');
    legend({'phi = 0','phi = 90 deg','sin(theta)'});
    title(sprintf('Hertzian dipole, D_{max} = %.3f', nf2ff.Dmax));
end

pass = ts_check_rel('directivity of a Hertzian dipole', nf2ff.Dmax, 1.5, limit_Dmax);

% the pattern, away from the nulls where a relative check says nothing
idx = (theta > 10/180*pi) & (theta < 170/180*pi);
for p=1:numel(phi)
    cut = abs(E_theta(idx,p)) / max(abs(E_theta(:,p)));
    pass = ts_check_rel(sprintf('sin(theta) pattern at phi = %g deg', phi(p)/pi*180), ...
                        cut, sin(theta(idx)), limit_pattern) && pass;
end

pass = ts_check_rel('pattern is rotationally symmetric', ...
                    abs(E_theta(idx,2)), abs(E_theta(idx,1)), limit_symmetry) && pass;

crosspol = max(abs(E_phi(:))) / max(abs(E_theta(:)));
pass = ts_check('no cross polarization', crosspol <= limit_crosspol, ...
                'max|E_phi| / max|E_theta| = %.3g (limit %.3g)', ...
                crosspol, limit_crosspol) && pass;

pass = ts_finish(opt, mfilename, pass, Sim_Path);
