function pass = Coax( openEMS_options, options )

physical_constants;


ENABLE_PLOTS = 1;
CLEANUP = 1;        % if enabled and result is PASS, remove simulation folder
STOP_IF_FAILED = 1; % if enabled and result is FAILED, stop with error
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

% LIMITS
upper_error = 0.03; % max +3%
lower_error = 0.01; % max -1%

% structure
length = 1000;
coax_rad_i  = 100;
coax_rad_ai = 230;
coax_rad_aa = 240;
mesh_res = [5 5 5];
f_start = 0;
f_stop = 1e9;

Sim_Path = 'tmp_Coax';
Sim_CSX = 'coax.xml';

[status,message,messageid]=rmdir(Sim_Path,'s');
[status,message,messageid]=mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(5000,1e-6);
FDTD = SetGaussExcite(FDTD,0,f_stop);
FDTD = SetBoundaryCond(FDTD,{'PEC','PEC','PEC','PEC','PEC','PML_8'});

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = -2.5*mesh_res(1)-coax_rad_aa : mesh_res(1) : coax_rad_aa+2.5*mesh_res(1);
mesh.y = mesh.x;
mesh.z = 0 : mesh_res(3) : length;
mesh.z = linspace(0,length,numel(mesh.z));
CSX = DefineRectGrid(CSX, 1e-3,mesh);

% create a perfect electric conductor
CSX = AddMetal(CSX,'PEC');

%%% coax
start = [0, 0 , 0];stop = [0, 0 , length];
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

% %dump
% CSX = AddDump(CSX,'Et_',0,2);
% start = [mesh.x(1) , 0 , mesh.z(1)];
% stop = [mesh.x(end) , 0 , mesh.z(end)];
% CSX = AddBox(CSX,'Et_',0 , start,stop);
% 
% CSX = AddDump(CSX,'Ht_',1,2);
% CSX = AddBox(CSX,'Ht_',0,start,stop);

%voltage calc
CSX = AddProbe(CSX,'ut1',0);
start = [ coax_rad_i 0 length/2 ];stop = [ coax_rad_ai 0 length/2 ];
CSX = AddBox(CSX,'ut1', 0 ,start,stop);

%current calc
CSX = AddProbe(CSX,'it1',1);
% mid = 0.5*(coax_rad_i+coax_rad_ai);
mid = coax_rad_i+3*mesh_res(1);
start = [ -mid -mid length/2 ];stop = [ mid mid length/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

% show structure
% CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

% run openEMS
folder = fileparts( mfilename('fullpath') );
Settings.LogFile = [folder '/' Sim_Path '/openEMS.log'];
Settings.Silent = SILENT;
RunOpenEMS( Sim_Path, Sim_CSX, openEMS_options, Settings );
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

% analytic formular for characteristic impedance
Z0 = sqrt(MUE0/EPS0) * log(coax_rad_ai/coax_rad_i) / (2*pi);
upper_limit = Z0 * (1+upper_error);
lower_limit = Z0 * (1-lower_error);

if ENABLE_PLOTS
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

pass = check_limits( Z, upper_limit, lower_limit );
if pass
    disp( 'combinedtests/Coax.m (characteristic impedance):  pass' );
else
    disp( 'combinedtests/Coax.m (characteristic impedance):  * FAILED *' );
end




if pass && CLEANUP
    rmdir( Sim_Path, 's' );
end
if ~pass && STOP_IF_FAILED
    error 'test failed';
end

