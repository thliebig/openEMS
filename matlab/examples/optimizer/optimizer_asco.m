%
% asco optimizer example -- optimize the turn number of a coil
%
% You need asco from http://asco.sf.net
% This is the main script.
%  - optimizer_simfun.m starts the simulator with a parameter set from
%    asco
%
% The goal is evaluated inside optimizer_simfun() to get as close to 2 uH.

% clear
clear
close all
clc

% setup the parameters
params = [];
params(end+1).name  = 'turns';
params(end).range   = [1 30];
params(end).value   = 4;
params(end).step    = 1;
params(end).active  = 1; % this parameter is to be optimized

% setup the simulation function
folder = fileparts( mfilename('fullpath') );
options.simfun = [folder '/optimizer_simfun.m'];

% additional options
% options.octave_exe = 'octave'; % must be newer than 3.2.4 (3.3.54 works)
options.clean = 1;

% start the optimization
[params_opt,result] = optimize( 'opttmp', params, options, 'asco' );

% display best value
disp( ['ASCO found the optimum turn number: ' num2str(params_opt(1).value) '  result: ' num2str(result)] );
