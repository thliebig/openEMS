%%%%%%%%%%%%%%%%%%%%%%%
% example demonstrating double drude meta-material
%
% tested with openEMS v0.0.28
%
% author: Thorsten Liebig @ 2010,2012
%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

%% setup the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
postproc_only = 0;      %set to 1 if the simulation is already done

Settings = [];
Settings.LogFile = 'openEMS.log';

pic_size = round([1400 1400/4]); %define the animation picture size

%simulation domain setup (in mm)
length = 500;
width = 10;
mesh_res = 0.5;        % mesh resolution
height = 3*mesh_res; % height is only 3 lines with PEC (top/bottom) --> quasi 2D

%FDTD setup
f0 = 5e9;         %center frequency
f_BW = f0/sqrt(2);  %bandwidth
MTM.eps_R = 1;
MTM.mue_R = 1;
MTM.f0 = f0;        %plasma frequency of the drude material
MTM.relaxTime = 5e-9; %relaxation time (smaller number results in greater losses, set to 0 to disable)
MTM.length = 250;  %length of the metamaterial
N_TS = 5e4;         %number of timesteps
endCriteria = 1e-5; %stop simulation if signal is at -50dB

%constants
physical_constants

%% define openEMS options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openEMS_opts = '-vvv';

Sim_Path = 'MTM_PW_Drude';
Sim_CSX = 'MTM_PW_Drude.xml';

if (postproc_only==0)
    
    if (exist(Sim_Path,'dir'))
        rmdir(Sim_Path,'s');
    end
    mkdir(Sim_Path);

    %% setup FDTD parameter & excitation function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FDTD = InitFDTD(N_TS,endCriteria,'OverSampling',10);
    FDTD = SetGaussExcite(FDTD,0,2*f0);
    BC = [1 1 0 0 2 2];
    FDTD = SetBoundaryCond(FDTD,BC);

    %% setup CSXCAD geometry & mesh %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CSX = InitCSX();
    mesh.x = -width/2 : mesh_res : width/2;
    mesh.y = -height/2 : mesh_res : height/2;
    mesh.z = -length/2 : mesh_res : length/2;
    CSX = DefineRectGrid(CSX, 1e-3,mesh);

    %% apply the plane wave excitation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    start=[-width/2 -height/2 ,mesh.z(3)];
    stop=[width/2 height/2 mesh.z(3)];
    CSX = AddExcitation(CSX,'excite',0,[0 1 0]); % excite E_y
    CSX = AddBox(CSX,'excite',0 ,start,stop);

    %% apply drude material %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CSX = AddLorentzMaterial(CSX,'drude');
    CSX = SetMaterialProperty(CSX,'drude','Epsilon',MTM.eps_R,'EpsilonPlasmaFrequency',MTM.f0,'EpsilonRelaxTime',MTM.relaxTime);
    CSX = SetMaterialProperty(CSX,'drude','Mue',MTM.mue_R,'MuePlasmaFrequency',MTM.f0,'MueRelaxTime',MTM.relaxTime);
    start=[mesh.x(1)   mesh.y(1)   -MTM.length/2];
    stop =[mesh.x(end) mesh.y(end)  MTM.length/2];
    CSX = AddBox(CSX,'drude', 10 ,start,stop);

    %% define dump boxes... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CSX = AddDump(CSX,'Et','FileType',1,'SubSampling','10,10,1');
    start = [mesh.x(2) ,0 , mesh.z(1)];
    stop = [mesh.x(end-1) , 0 , mesh.z(end)];
    CSX = AddBox(CSX,'Et',0 , start,stop);

    %% Write openEMS compatoble xml-file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    WriteOpenEMS([Sim_Path '/' Sim_CSX],FDTD,CSX);

    %% run openEMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RunOpenEMS(Sim_Path, Sim_CSX, openEMS_opts, Settings);

end

%% plot the drude type material dependency
f = linspace(0.1*f0,2*f0,501);
w = 2*pi*f;
epsr = MTM.eps_R * (1 - (2*pi*MTM.f0)^2./( w.^2 - 1j*w./MTM.relaxTime ));
muer = MTM.mue_R * (1 - (2*pi*MTM.f0)^2./( w.^2 - 1j*w./MTM.relaxTime ));
plot(f,real(epsr),'Linewidth',2);
hold on
grid on
plot(f,imag(epsr),'r--','Linewidth',2);
plot(f,real(muer),'c-.','Linewidth',2);
plot(f,imag(muer),'m-.','Linewidth',2);
ylim([-10 MTM.eps_R])
% l=legend('\Re \epsilon_r','\Im \epsilon_r','\Re \mue_r','\Im \mue_r');
l=legend('$\Re\{\varepsilon_r\}$','$\Im\{\varepsilon_r\}$','$\Re\{\mu_r\}$','$\Im\{\mu_r\}$');
set(l,'Interpreter','latex','Fontsize',12)

%% plot E-fields
freq = [f0/sqrt(2) f0 f0*sqrt(2)];
field = ReadHDF5FieldData([Sim_Path '/Et.h5']);
mesh_h5 = ReadHDF5Mesh([Sim_Path '/Et.h5']);

ET = ReadUI('et',Sim_Path);
ef = DFT_time2freq(ET.TD{1}.t,ET.TD{1}.val,freq);

field_FD = GetField_TD2FD(field, freq);

mesh.x = linspace(-500,500,numel(mesh_h5.lines{1})); %make animation wider...
mesh.y = mesh_h5.lines{2};
mesh.z = mesh_h5.lines{3};

[X Z] = meshgrid(mesh.x,mesh.z);
X = X';
Z = Z';

for n=1:numel(field_FD.FD.values)
    Ec{n} = squeeze(field_FD.FD.values{n}/ef(n));
end

%%
figure('Position',[10 100 pic_size(1) pic_size(2)]);
phase = linspace(0,2*pi,21);
disp('press CTRL+C to stop animation');
while (1)
    for ph = phase(1:end-1)
        for n=1:numel(Ec)
            subplot(1,numel(Ec),n)
            E = real(Ec{n}.*exp(1j*ph));
            surf(X,Z,E(:,:,2));
            title(['f_0 = ' num2str(freq(n)*1e-9) ' GHz'])
        end
        pause(0.1);
    end
end
