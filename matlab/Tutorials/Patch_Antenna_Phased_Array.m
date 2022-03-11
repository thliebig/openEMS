%
% Tutorials / Patch Antenna Phased Array
%
% Description at:
%
% Tested with
%  - Matlab 2011a
%  - Octave 4.0
%  - openEMS v0.0.33
%
% References:
% [1] Y. Yusuf and X. Gong, “A low-cost patch antenna phased array with
%   analog beam steering using mutual coupling and reactive loading,” IEEE
%   Antennas Wireless Propag. Lett., vol. 7, pp. 81–84, 2008.
% [2] S. Otto, S. Held, A. Rennings, and K. Solbach,
%   "Array and multiport antenna farfield simulation using
%   EMPIRE, MATLAB and ADS," 39th European Microwave Conf. (EuMC 2009),
%   Sept. 29 – Oct. 1, Rome, Italy, pp. 1547-1550, 2009.
% [3] K. Karlsson, J. Carlsson, I. Belov, G. Nilsson, and P.-S. Kildal,
%   “Optimization of antenna diversity gain by combining full-wave and
%   circuit simulations,” in Proc. Second European Conference on Antennas
%   and Propagation EuCAP 2007, 11–16 Nov. 2007, pp. 1–5.
%
% (C) 2013-2015 Thorsten Liebig <thorsten.liebig@gmx.de>


close all
clear
clc

% we need the "Cuircuit Toolbox"
addpath('C:\CTB');
% get the latest version from:
% using git: https://github.com/thliebig/CTB
% or zip: https://github.com/thliebig/CTB/archive/master.zip

% set this to 0 to NOT run a reference simulation with the given C2 and C3
% for comparison
do_reference_simulation = 1;

% set to 1 if you want to run AppCSXCAD to see the simulated structure
show_structure = 1;

% set this to 1, to force openEMS to run again even if the data already exist
force_rerun = 0;

% frequency range of interest
f = linspace( 1e9, 5e9, 1601 );

% resonant frequency for far-field calculations
f0 = 3e9;

% capacities for port 2 and 3 to shift the far-field pattern
C2 = 0.2e-12;
C3 = 0.2e-12;

Sim_Path_Root = ['tmp_' mfilename];

%% calculate the full S-parameter set for all 3 patch antennas running 3
% individual openEMS simulations in which one antenna is active and the
% other two a passive (50 Ohm load) respectively
xpos = [0 -41 41]; % x-center position of the 3 antennas
caps = [0 0 0];
resist = [50 50 50];

spara = [];
color_code = {'k-','r--','m-.'};

for n=1:3
    active = [0 0 0];
    active(n) = 1; % activate antenna n
    Sim_Path = [Sim_Path_Root '_' num2str(n)]; % create an individual path
    [port{n} nf2ff{n}] = Patch_Antenna_Array(Sim_Path, ((exist(Sim_Path,'dir')>0) && (force_rerun==0)), show_structure, xpos, caps, resist, active);
    port{n} = calcPort( port{n}, Sim_Path, f, 'RefImpedance', 50);
    nf2ff{n} = CalcNF2FF(nf2ff{n}, Sim_Path, f0, [-180:2:180]*pi/180, 0);

    figure
    hold on
    grid on
    for p=1:3
        I(p,n) = interp1(f, port{n}{p}.if.tot,f0);
        P_in(p) = 0.5*interp1(f, port{n}{n}.uf.inc,f0)*conj(interp1(f, port{n}{n}.if.inc,f0));
        spara(p,n,:) = port{n}{p}.uf.ref./ port{n}{n}.uf.inc;
        plot(f, squeeze(20*log10(abs(spara(p,n,:)))),color_code{p},'Linewidth',2);
    end
end

%% export sparameter to touchstone file
write_touchstone('s',f,spara,[Sim_Path_Root '.s3p']);

% instructions for Qucs:
% load the written touchstone file
% attach C2 and C3 to port 2 and 3
% attach a signal port to port 1
% probe the currents going into port 1 to 3

% example currents for ports 1 to 3 for C2 = 0.2pF and C3=0.2pF
I_qucs(1,1) = 0.00398-0.000465j;
I_qucs(2,1) = 2.92e-5-0.000914j;
I_qucs(3,1) = 2.92e-5-0.000914j;

disp(['I2/I1: Qucs: ' num2str(I_qucs(2)/I_qucs(1)) ' (defined manually)'])
disp(['I3/I1: Qucs: ' num2str(I_qucs(3)/I_qucs(1)) ' (defined manually)'])

%% Calculate the currents of port 1 to 3 using Matlab [1]
z = s2z(spara);

Z2 = 1/(1j*2*pi*f0*C2);
Z3 = 1/(1j*2*pi*f0*C3);

z23(1,1) = interp1(f,squeeze(z(2,2,:)),f0) + Z2;
z23(1,2) = interp1(f,squeeze(z(2,3,:)),f0);
z23(2,1) = interp1(f,squeeze(z(3,2,:)),f0);
z23(2,2) = interp1(f,squeeze(z(3,3,:)),f0) + Z3;

%set input/feeding current of port 1 to 1mA
I_out(1,1) = 1e-3;
% calc current for port 2 and 3
I_out([2 3],1) = z23\[-interp1(f,squeeze(z(2,1,:)),f0);-interp1(f,squeeze(z(3,1,:)),f0)]*I_out(1);

disp(['I2/I1: Matlab: ' num2str(I_out(2)/I_out(1))])
disp(['I3/I1: Matlab: ' num2str(I_out(3)/I_out(1))])


%% do a reference simulation for the given C2/C3 values
if (do_reference_simulation)
    active = [1 0 0];
    caps = [0 C2 C3];
    resist = [50 inf inf];
    Sim_Path = [Sim_Path_Root '_C2=' num2str(C2*1e12) '_C3=' num2str(C3*1e12)];
    [port_ref nf2ff_ref] = Patch_Antenna_Array(Sim_Path, ((exist(Sim_Path,'dir')>0) && (force_rerun==0)), show_structure, xpos, caps, resist, active);
    port_ref = calcPort( port_ref, Sim_Path, f, 'RefImpedance', 50);
    nf2ff_ref = CalcNF2FF(nf2ff_ref, Sim_Path, f0, [-180:2:180]*pi/180, 0);

    % extract currents from reference simulation
    for p=1:3
        I_ref(p,1) = interp1(f, port_ref{p}.if.tot,f0);
    end

    disp(['I2/I1: openEMS: ' num2str(I_ref(2)/I_ref(1))])
    disp(['I3/I1: openEMS: ' num2str(I_ref(3)/I_ref(1))])
end

%% calculate and apply weighting coefficients [3]
% calculate
coeff = I\I_out;

% apply
E_ff_phi = 0*nf2ff{1}.E_phi{1};
E_ff_theta = 0*nf2ff{1}.E_phi{1};
for n=1:3
    E_ff_phi = E_ff_phi + coeff(n)*nf2ff{n}.E_phi{1};
    E_ff_theta = E_ff_theta + coeff(n)*nf2ff{n}.E_theta{1};
end

%% plot far-field patterns
figure
polar([-180:2:180]'*pi/180,abs(E_ff_phi(:))/max(abs(E_ff_phi(:))));
hold on
if (do_reference_simulation)
    polar([-180:2:180]'*pi/180,abs(nf2ff_ref.E_norm{1}(:,1))/max(abs(nf2ff_ref.E_norm{1}(:,1))),'r--');
end
title('normalized far-field pattern','Interpreter', 'none')
legend('calculated','reference')


