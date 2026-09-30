function FDTD = SetDiracExcite(FDTD)
% function FDTD = SetDiracExcite(FDTD)
%
% Set a dirac pulse as excitation signal: the full amplitude is applied
% for one timestep and zero in the next.
%
% The pulse is two timesteps long, so its duration -- and with it its
% bandwidth -- depends on the mesh, and it is not band-limited. Unlike
% SetGaussExcite this function sets no maximum frequency for the
% simulation, so probes and field dumps are sampled every timestep: with
% field dumps enabled that is slow and writes very large files, so limit
% 'NrTS' in InitFDTD first. Prefer SetGaussExcite unless an impulse
% response is really what is wanted.
%
% see also SetGaussExcite SetSinusExcite SetStepExcite SetCustomExcite
%
% e.g
%
%     FDTD = InitFDTD('NrTS', 1000);
%     FDTD = SetDiracExcite(FDTD);
%
% openEMS matlab interface
% -----------------------
% author: Thorsten Liebig

FDTD.Excitation.ATTRIBUTE.Type=2;
