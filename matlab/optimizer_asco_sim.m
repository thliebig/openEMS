function optimizer_asco_sim( optimdir, inputfile, outputfile, simfun )
%optimizer_asco_sim( optimdir, inputfile, outputfile, simfun )
%
% This function is called from general.sh. Do not call it yourself.
%
% tasks:
%  - set correct matlab path
%  - evaluate input file
%  - start simulation or get result from cache
%  - post-process simulation results
%  - create output file (important: needs single \n at the first line and double \n at the last line!)

error( nargchk(4,4,nargin) );

% add CSXCAD and openEMS to the matlab path
folder = fileparts( mfilename('fullpath') );
addpath( folder );
addpath( [folder '/../../CSXCAD/matlab'] );

% change to optimdir
olddir = pwd;
cd( optimdir );

% read parameters set by asco
if ~isempty( strfind(inputfile,'-') )
    % matlab cannot execute a file with dashes...
    inputfile2 = strrep( inputfile,'-','_' );
    movefile( [inputfile '.m'], [inputfile2 '.m'] );
    run( inputfile2 );
    movefile( [inputfile2 '.m'], [inputfile '.m'] );
end
% now a structure named 'params' is available

% check cache
folder = create_folder_name( params );
if exist( ['./' folder], 'dir' ) && exist( ['./' folder '/result.mat'], 'file' )
    % read cache
    disp( 'CACHE HIT' );
    result = load( [folder '/result.mat'], 'result' );
    result = result.result;
else
    % start simulation in folder <folder>
    disp( ['starting simulation function ' simfun] );
    disp( ['         simulation folder   ' folder] );
    [simfun_folder,simfun] = fileparts(simfun);
    oldpath = path;
    addpath( simfun_folder );
    fhandle = str2func(simfun); % does not work for octave-3.2.4!
    path( oldpath );
    mkdir( folder );
    result = fhandle(folder,params);
    save( [folder '/result.mat'], 'result', '-mat' );
end

% write results for asco
fid = fopen( outputfile, 'wt' );
fprintf( fid, '\nvalue= %e\n\n', result );
fclose( fid );

% update best result
best = [];
best.result = result;
best.params = params;
if exist( [pwd '/best_result.mat'], 'file' )
    old = load( 'best_result.mat', 'best' );
    if old.best.result > best.result
        save( 'best_result.mat', 'best', '-mat' );
    end
else
    save( 'best_result.mat', 'best', '-mat' );
end

% restore old folder
cd( olddir );







function folder = create_folder_name( params )
params = orderfields( params );
folder = 'opt';
fnames = fieldnames(params);
for n=1:numel(fnames)
    folder = [folder '_' fnames{n} '=' num2str(params.(fnames{n}))];
end
