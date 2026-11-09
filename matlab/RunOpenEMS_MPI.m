function RunOpenEMS_MPI(Sim_Path, Sim_File, opts, Settings)
% function RunOpenEMS_MPI(Sim_Path, Sim_File, opts, Settings)
%
% Run an openEMS simulation on Settings.MPI.NrProc MPI ranks, see RunOpenEMS
%
%     Settings.MPI.NrProc = 4;
%     Settings.MPI.Binary = '/opt/openEMS/bin/openEMS';   % default: openEMS on the PATH
%     Settings.MPI.Launcher = 'mpiexec';                  % default, or e.g. 'srun'
%     Settings.MPI.GlobalArgs = '';                       % arguments of the launcher
%     Settings.MPI.Hosts = {'host1','host2'};             % optional
%     RunOpenEMS(Sim_Path, Sim_File, '--mpi-split=2,2,1', Settings)
%
% See also RunOpenEMS
%
% openEMS matlab interface
% -----------------------

if ~isunix
    error('openEMS:RunOpenEMS_MPI','MPI runs are available on unix only');
end

mpi = Settings.MPI;
if ~isfield(mpi,'Binary')
    mpi.Binary = 'openEMS';
end
if ~isfield(mpi,'Launcher')
    mpi.Launcher = 'mpiexec';
end
if ~isfield(mpi,'GlobalArgs')
    mpi.GlobalArgs = '';
end
hosts = '';
if isfield(mpi,'Hosts')
    hosts = [' -host ' strjoin(mpi.Hosts, ',')];
end

command = [mpi.Launcher ' ' mpi.GlobalArgs hosts ' -n ' int2str(mpi.NrProc) ' "' mpi.Binary '" ' Sim_File ' ' opts];
if isfield(Settings,'LogFile')
    command = [command ' > "' Settings.LogFile '" 2>&1'];
end

savePath = pwd;
cd(Sim_Path);
status = system(command);
cd(savePath);
if (status~=0)
    error('openEMS:RunOpenEMS_MPI', ['the MPI run of openEMS failed, exit status ' num2str(status)]);
end
