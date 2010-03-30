function pass = cavity

physical_constants;


ENABLE_PLOTS = 1;
CLEANUP = 0;        % if enabled and result is PASS, remove simulation folder
STOP_IF_FAILED = 1; % if enabled and result is FAILED, stop with error

% LIMITS
upper_error = 0.036; % max +3.6%
lower_error = 0;     % max -0.0%

% structure
a = 5e-2;
b = 2e-2;
d = 6e-2;
if ~((b<a) && (a<d))
    error 'correct the dimensions of the cavity'
end

f_start = 0;
f_stop = 10e9;

Sim_Path = 'tmp';
Sim_CSX = 'cavity.xml';

[status,message,messageid]=mkdir(Sim_Path);

%setup FDTD parameter
FDTD = InitFDTD( 10000,1e-6 );
FDTD = SetGaussExcite(FDTD,(f_stop-f_start)/2,(f_stop-f_start)/2);
BC = [0 0 0 0 0 0]; % PEC boundaries
FDTD = SetBoundaryCond(FDTD,BC);

%setup CSXCAD geometry
CSX = InitCSX();
mesh.x = linspace(0,a,35);
mesh.y = linspace(0,b,35);
mesh.z = linspace(0,d,35);
CSX = DefineRectGrid(CSX, 1,mesh);

% excitation
pos1 = [mesh.x(10) mesh.y(10) mesh.z(10)];
pos2 = [mesh.x(10) mesh.y(11) mesh.z(10)];
CSX = AddExcitation(CSX,'excite',1,[0 1 0]);
CSX = AddBox(CSX, 'excite', 1, pos1, pos2);
 
% %dump
% CSX = AddDump(CSX,'Et_',0,2);
% pos1 = [mesh.x(1) mesh.y(10) mesh.z(1)];
% pos2 = [mesh.x(end) mesh.y(10) mesh.z(end)];
% CSX = AddBox(CSX,'Et_',0 , pos1,pos2);
% 
% %dump
% CSX = AddDump(CSX,'Et2_',0,2);
% pos1 = [mesh.x(1) mesh.y(1) mesh.z(1)];
% pos2 = [mesh.x(end) mesh.y(1) mesh.z(end)];
% CSX = AddBox(CSX,'Et2_',0 , pos1,pos2);
% 
% %dump
% CSX = AddDump(CSX,'Et3_',0,2);
% pos1 = [mesh.x(1) mesh.y(end-1) mesh.z(1)];
% pos2 = [mesh.x(end) mesh.y(end-1) mesh.z(end)];
% CSX = AddBox(CSX,'Et3_',0 , pos1,pos2);

%voltage calc
pos1 = [mesh.x(15) mesh.y(15) mesh.z(15)];
pos2 = [mesh.x(15) mesh.y(16) mesh.z(15)];
CSX = AddProbe(CSX,'ut1y',0);
CSX = AddBox(CSX,'ut1y', 0 ,pos1,pos2);
pos1 = [mesh.x(floor(end/2)) mesh.y(floor(end/2)) mesh.z(end)];
pos2 = [mesh.x(floor(end/2)) mesh.y(floor(end/2)) mesh.z(end-1)];
CSX = AddProbe(CSX,'ut1z',0);
CSX = AddBox(CSX,'ut1z', 0 ,pos1,pos2);

%Write openEMS compatible xml-file
WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

%cd to working dir and run openEMS
savePath = pwd();
cd(Sim_Path); %cd to working dir
invoke_openEMS( Sim_CSX );
UI = ReadUI( {'ut1x', 'ut1y', 'ut1z'} );
cd(savePath);



%
% analysis
%

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
f_TE103 = c0/(2*pi) * k(1,0,3);
f_TE201 = c0/(2*pi) * k(2,0,1);
f_TE202 = c0/(2*pi) * k(2,0,2);
f_TM110 = c0/(2*pi) * k(1,1,0);
f_TM210 = c0/(2*pi) * k(2,1,0);

f_TE = [f_TE101 f_TE102 f_TE103 f_TE201 f_TE202];
f_TM = [f_TM110 f_TM210];

if ENABLE_PLOTS
    figure
    plot(f/1e9,abs(uy))
    max1 = max(abs(uy));
    hold on
    plot( repmat(f_TE,2,1)/1e9, repmat([0; max1],1,numel(f_TE)), 'm-.', 'LineWidth', 2 )
    xlabel('Frequency (GHz)')
    legend( {'u_y','theoretical'} );

    figure
    plot(f/1e9,abs(uz))
    max1 = max(abs(uz));
    hold on
    plot( repmat(f_TM,2,1)/1e9, repmat([0; max1],1,numel(f_TM)), 'm-.', 'LineWidth', 2 )
    xlabel('Frequency (GHz)')
    legend( {'u_z','theoretical'} );
end

pass = 1;
% pass = check_limits( Z, upper_limit, lower_limit );
% if pass
%     disp( 'combinedtests/Coax.m (characteristic impedance):  pass' );
% else
%     disp( 'combinedtests/Coax.m (characteristic impedance):  * FAILED *' );
% end




if pass && CLEANUP
    rmdir( [Sim_Path '/' Sim_CSX], 's' );
end
if ~pass && STOP_IF_FAILED
    error 'test failed';
end
