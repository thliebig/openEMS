function OptimizeCondSheetParameter
% function OptimizeCondSheetParameter
% 
% internal openEMS function to create the "cond_sheet_parameter.h" header
% file containing optimized parameter values for the conducting sheet model
%
% (c) 2012: Thorsten Liebig, thorsten.liebig@gmx.de


close all
clear
clc


X = [1 1 1 1 1];
lb = zeros(size(X));
ub = ones(size(X))*1000;


Omega = linspace(0,20,51);
Omega = Omega(2:end);

Omega_fine = linspace(0,20,15001);
Omega_fine = Omega_fine(2:end);

options = optimset('fzero');
% options = optimset(options,'Display','iter')
options = optimset(options,'MaxFunEvals',5000);
options = optimset(options,'MaxIter',5000);

omg_stop_str = [];
omg_critical_20 = [];
omg_critical_5 = [];
g_str = [];
r1_str= [];
r2_str= [];
l1_str= [];
l2_str= [];

numPara = 30;
factor = 1.5;

for p = 1:numPara  
    
    X = lsqnonlin(@(X)CalcYDiff(Omega,X),X,lb,ub,options);

    omg_stop_str = [omg_stop_str  num2str(Omega(end),10) ','];
    g_str = [g_str num2str(X(1),10) ','];
    r1_str = [r1_str num2str(X(2),10) ','];
    l1_str = [l1_str num2str(X(3),10) ','];
    r2_str = [r2_str num2str(X(4),10) ','];
    l2_str = [l2_str num2str(X(5),10) ','];

    Ys = tanh((1+1j)*sqrt(Omega_fine))/(1+1j)./sqrt(Omega_fine);

    err = (CalcYDiff(Omega_fine,X))./real(1./Ys);
    
    fc = Omega_fine(find(abs(err)>0.2,1,'last'));
    if (isempty(fc))
        fc = 0;
    end
    omg_critical_20 = [omg_critical_20 num2str(fc,10) ','];
    
    fc = Omega_fine(find(abs(err)>0.05,1,'last'));
    if (isempty(fc))
        fc = 0;
    end
    omg_critical_5 = [omg_critical_5 num2str(fc,10) ','];
    
%     disp(['max error: ' num2str(max(abs(err)*100)) ])
%     figure
%     plot(Omega_fine,err*100,'g--');
    
    Omega_fine = Omega_fine*factor;
    Omega = Omega*factor;
end

%% write to file
fid = fopen('cond_sheet_parameter.h','w');

fprintf(fid,'// This is a list of conducting sheet model parameter for different ranges of omega = w/w0\n');
fprintf(fid,'// This file was created automatically using Matlab: OptimizeCondSheetParameter.m \n');
fprintf(fid,'// Do not change this file! \n');
fprintf(fid,'// Creation: %s \n\n',datestr(now));



fprintf(fid,'unsigned int numOptPara=%d;\n',numPara);

fprintf(fid,'double omega_stop[%d]={%s};\n',numPara,omg_stop_str(1:end-1));
fprintf(fid,'double omega_critical_5[%d]={%s};\n',numPara,omg_critical_5(1:end-1));
fprintf(fid,'double omega_critical_20[%d]={%s};\n',numPara,omg_critical_20(1:end-1));
fprintf(fid,'double g[%d]={%s};\n',numPara,g_str(1:end-1));
fprintf(fid,'double r1[%d]={%s};\n',numPara,r1_str(1:end-1));
fprintf(fid,'double l1[%d]={%s};\n',numPara,l1_str(1:end-1));
fprintf(fid,'double r2[%d]={%s};\n',numPara,r2_str(1:end-1));
fprintf(fid,'double l2[%d]={%s};\n',numPara,l2_str(1:end-1));

fclose(fid);

function Ydiff =  CalcYDiff(omega, X)

Ys = tanh((1+1j)*sqrt(omega))/(1+1j)./sqrt(omega);

Y = X(1) + 1./(X(2)+1j*omega*X(3)) + 1./(X(4)+1j*omega*X(5));

Ydiff = real(1./Ys)-real(1./Y);