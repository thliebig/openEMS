function pass = Coax(varargin)
%pass = Coax(<key>, <value>, ...)
%
% Coaxial line in cartesian coordinates, excited with the analytic TEM mode
% weighting. Checks the characteristic impedance against
%
%     Z0 = sqrt(mue0/eps0) * log(r_outer/r_inner) / (2*pi)
%
% The limits are asymmetric: the staircased round conductors raise the
% simulated impedance slightly above the analytic value of the ideal line.
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
upper_error = 0.03; % max +3%
lower_error = 0.01; % max -1%

% structure
coax_length = 1000;
coax_rad_i  = 100;
coax_rad_ai = 230;
coax_rad_aa = 240;
mesh_res = [5 5 5];
f_start = 0;
f_stop = 1e9;

Sim_Path = ts_sim_path(mfilename('fullpath'));
Sim_CSX = 'coax.xml';

%setup FDTD parameter
% -40 dB is the floor this geometry reaches: the PEC short at z=0 turns the
% residual DC of the excitation into a standing current the PML cannot absorb.
% 1e-4 is where Z has converged (50.49..51.12 Ohm, same as 1e-5 and 1e-6);
% anything stricter only runs into NrTS.
FDTD = InitFDTD('NrTS', 5000, 'EndCriteria', 1e-4);
FDTD = SetGaussExcite(FDTD,f_stop/2,f_stop/3);
FDTD = SetBoundaryCond(FDTD,{'PEC','PEC','PEC','PEC','PEC','PML_8'});

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = -2.5*mesh_res(1)-coax_rad_aa : mesh_res(1) : coax_rad_aa+2.5*mesh_res(1);
mesh.y = mesh.x;
mesh.z = 0 : mesh_res(3) : coax_length;
mesh.z = linspace(0,coax_length,numel(mesh.z));
CSX = DefineRectGrid(CSX, 1e-3,mesh);

% create a perfect electric conductor
CSX = AddMetal(CSX,'PEC');

%%% coax
start = [0, 0 , 0];stop = [0, 0 , coax_length];
CSX = AddCylinder(CSX,'PEC',1 ,start,stop,coax_rad_i); % inner conductor
CSX = AddCylindricalShell(CSX,'PEC',0 ,start,stop,0.5*(coax_rad_aa+coax_rad_ai),(coax_rad_aa-coax_rad_ai)); % outer conductor

%%% add excitation
start(3) = 0; stop(3)=mesh_res(1)/2;
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = '(x)/(x*x+y*y)';
weight{2} = 'y/pow(rho,2)';
weight{3} = '0';
CSX = SetExcitationWeight(CSX, 'excite', weight );
CSX = AddCylindricalShell(CSX,'excite',0 ,start,stop,0.5*(coax_rad_i+coax_rad_ai),(coax_rad_ai-coax_rad_i));

%voltage calc
CSX = AddProbe(CSX,'ut1',0);
start = [ coax_rad_i 0 coax_length/2 ];stop = [ coax_rad_ai 0 coax_length/2 ];
CSX = AddBox(CSX,'ut1', 0 ,start,stop);

%current calc
CSX = AddProbe(CSX,'it1',1);
mid = coax_rad_i+3*mesh_res(1);
start = [ -mid -mid coax_length/2 ];stop = [ mid mid coax_length/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

% run openEMS
Settings.LogFile = [Sim_Path '/openEMS.log'];
Settings.Silent = opt.Silent;
RunOpenEMS( Sim_Path, Sim_CSX, opt.openEMS_opts, Settings );
UI = ReadUI( {[Sim_Path '/ut1'], [Sim_Path '/it1']} );


%
% analysis
%

f = UI.FD{2}.f;
u = UI.FD{1}.val;
i = UI.FD{2}.val;

f_idx_start = interp1( f, 1:numel(f), f_start, 'nearest' );
f_idx_stop  = interp1( f, 1:numel(f), f_stop,  'nearest' );
f = f(f_idx_start:f_idx_stop);
u = u(f_idx_start:f_idx_stop);
i = i(f_idx_start:f_idx_stop);

Z = abs(u./i);

% analytic formula for the characteristic impedance
Z0 = sqrt(MUE0/EPS0) * log(coax_rad_ai/coax_rad_i) / (2*pi);
upper_limit = Z0 * (1+upper_error);
lower_limit = Z0 * (1-lower_error);

if opt.Plots
    upper = upper_limit * ones(1,size(Z,2));
    lower = lower_limit * ones(1,size(Z,2));
    Z0_plot = Z0 * ones(1,size(Z,2));
    figure
    plot(f/1e9,[Z;upper;lower])
    hold on
    plot(f/1e9,Z0_plot,'m-.','LineWidth',2)
    hold off
    xlabel('Frequency (GHz)')
    ylabel('Impedance (Ohm)')
    legend( {'sim', 'upper limit', 'lower limit', 'theoretical'} );
end

pass = ts_check('characteristic impedance', ...
                check_limits( Z, upper_limit, lower_limit ), ...
                'Z = %.2f .. %.2f Ohm, limits %.2f .. %.2f Ohm (Z0 = %.2f Ohm)', ...
                min(Z), max(Z), lower_limit, upper_limit, Z0);

pass = ts_finish(opt, mfilename, pass, Sim_Path);
