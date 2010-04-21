function pass = Coax

physical_constants;


ENABLE_PLOTS = 1;
CLEANUP = 1;        % if enabled and result is PASS, remove simulation folder
STOP_IF_FAILED = 1; % if enabled and result is FAILED, stop with error

% LIMITS
upper_error = 0.036; % max +3.6%
lower_error = 0;     % max -0.0%

% structure
abs_length = 250;
length = 1000;
coax_rad_i  = 100;
coax_rad_ai = 230;
coax_rad_aa = 240;
mesh_res = [5 5 5];
f_start = 0;
f_stop = 1e9;

Sim_Path = 'tmp';
Sim_CSX = 'coax.xml';

[status,message,messageid]=rmdir(Sim_Path,'s');
[status,message,messageid]=mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD(5e5,1e-6);
FDTD = SetGaussExcite(FDTD,(f_stop-f_start)/2,(f_stop-f_start)/2);
BC = [1 1 1 1 1 1] * 0;
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = -2.5*mesh_res(1)-coax_rad_aa : mesh_res(1) : coax_rad_aa+2.5*mesh_res(1);
mesh.y = mesh.x;
mesh.z = 0 : mesh_res(3) : length;
mesh.z = linspace(0,length,numel(mesh.z) + 4-mod(numel(mesh.z),4)); % make it compatible with sse-engine
CSX = DefineRectGrid(CSX, 1e-3,mesh);

%create copper
CSX = AddMetal(CSX,'PEC');

%%%fake pml
finalKappa = 0.3/abs_length^4;
finalSigma = finalKappa*MUE0/EPS0;
CSX = AddMaterial(CSX,'pml');
CSX = SetMaterialProperty(CSX,'pml','Kappa',finalKappa);
CSX = SetMaterialProperty(CSX,'pml','Sigma',finalSigma);
CSX = SetMaterialWeight(CSX,'pml','Kappa',['pow(abs(z)-' num2str(length-abs_length) ',4)']);
CSX = SetMaterialWeight(CSX,'pml','Sigma',['pow(abs(z)-' num2str(length-abs_length) ',4)']);

%%% coax
start = [0, 0 , 0];stop = [0, 0 , length];
CSX = AddCylinder(CSX,'PEC',0 ,start,stop,coax_rad_i); % inner conductor
CSX = AddCylindricalShell(CSX,'PEC',0 ,start,stop,0.5*(coax_rad_aa+coax_rad_ai),(coax_rad_aa-coax_rad_ai)); % outer conductor

%%% add PML
start(3) = length-abs_length;
CSX = AddCylindricalShell(CSX,'pml',0 ,start,stop,0.5*(coax_rad_i+coax_rad_ai),(coax_rad_ai-coax_rad_i));
start(3) = 0; stop(3)=mesh_res(1)/2;
CSX = AddExcitation(CSX,'excite',0,[1 1 0]);
weight{1} = '(x)/(x*x+y*y)';
weight{2} = 'y/pow(rho,2)';
weight{3} = 0;
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
mid = 0.5*(coax_rad_i+coax_rad_ai);
start = [ -mid -mid length/2 ];stop = [ mid mid length/2 ];
CSX = AddBox(CSX,'it1', 0 ,start,stop);

%Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
invoke_openEMS( Sim_CSX );
UI = ReadUI( {'ut1','it1'} );
cd(savePath);



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
