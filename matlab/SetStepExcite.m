function FDTD = SetStepExcite(FDTD)
% function FDTD = SetStepExcite(FDTD)
%
% Set a heaviside step as excitation signal: the full amplitude is
% applied from the first timestep on and never returns to zero.
%
% Because the excitation never ends, the energy in the simulation domain
% does not decay and the end criteria is never reached: limit the run
% with 'NrTS' in InitFDTD or with MaxTime. The step is not band-limited
% either, and as for SetDiracExcite no maximum frequency is set for the
% simulation, so probes and field dumps are sampled every timestep. For a
% step with a finite rise time use SetCustomExcite.
%
% see also SetGaussExcite SetSinusExcite SetDiracExcite SetCustomExcite
%
% e.g
%
%     FDTD = InitFDTD('NrTS', 10000);
%     FDTD = SetStepExcite(FDTD);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.Excitation.ATTRIBUTE.Type=3;
