function [params,result] = optimize( optimdir, params, options, algorithm )
%params = optimize( optimdir, params, options, algorithm )
%
% input:
%   optimdir:   folder where to optimize
%   params:     array of structures with parameters to optimize
%     .name       char string   parameter name (e.g. 'length1')
%     .value      number  
%     .step       number        discretization of .value
%     .range      row vector    range of .value (e.g. [1 10])
%     .active     (0/1)         1=optimize this parameter
%   options:    structure
%     .folder_matlabstart   set the startup folder for matlab or octave
%     .simfun               char string with the simulation function name
%     .octave_exe           if this field is present, octave is used
%     .clean                clean the optimization folder before optimization start
%   algorithm:  'asco' or 'simplex-downhill'
%   
% output:
%   params: optimal values
%   result: optimization criterion for optimal values
%
% example:
%   see openEMS/matlab/examples/optimizer
%
% notes:
%   Create a file named 'STOP_OPTIMIZATION' in the optimdir folder to stop
%   the optimization process (or press Ctrl+C).
%
% (C) 2010 Sebastian Held <sebastian.held@gmx.de>

error( nargchk(3,4,nargin) );

% default to simplex-downhill
if nargin < 4
    algorithm = 'simplex-downhill';
end
if ~strcmp( algorithm, 'asco' )
    algorithm = 'simplex-downhill';
end

optimdir = absolutepath( optimdir );
if isfield(options,'clean') && (options.clean == 1)
    [a,a,a] = rmdir( optimdir, 's' );
end
[a,a,a] = mkdir( optimdir );
oldfolder = cd( optimdir );

if strcmp( algorithm, 'asco' )
    % ---------------------------------------------------------------------
    % selected algorithm: ASCO
    % http://asco.sourceforge.net/

%     if ~exist( 'asco', 'file' ) && ~exist( 'asco.exe', 'file' )
%         error 'asco was not found in PATH. Download from http://asco.sf.net/'
%     end
    
    % create asco config file
    fid = fopen( 'asco.cfg', 'wt' );
    fprintf( fid, '* asco configuration file\n' );
    fprintf( fid, '* http://asco.sourceforge.net/\n' );
    fprintf( fid, '* \n' );
    fprintf( fid, '* Optimization for openEMS\n\n' );
    fprintf( fid, '#Optimization Flow#\n' );
    fprintf( fid, 'Alter:no            $do we want to do corner analysis?\n' );
    fprintf( fid, 'MonteCarlo:no       $do we want to do MonteCarlo analysis?\n' );
    fprintf( fid, 'AlterMC cost:0.00   $point below which ALTER and/or MONTECARLO can start\n' );
    fprintf( fid, 'ExecuteRF:no        $Execute or no the RF module to add RF parasitics?\n' );
    fprintf( fid, '#\n\n' );
    fprintf( fid, '#DE#\n' );
    fprintf( fid, 'choice of method:3\n' );
    fprintf( fid, 'maximum no. of iterations:50\n' );
    fprintf( fid, 'Output refresh cycle:2\n' );
    fprintf( fid, 'No. of parents NP:10\n' );
    fprintf( fid, 'Constant F:0.85\n' );
    fprintf( fid, 'Crossing Over factor CR:1\n' );
    fprintf( fid, 'Seed for pseudo random number generator:3\n' );
    fprintf( fid, 'Minimum Cost Variance:1e-6\n' );
    fprintf( fid, 'Cost objectives:10\n' );
    fprintf( fid, 'Cost constraints:100\n' );
    fprintf( fid, '#\n\n' );
    
    fprintf( fid, '# Parameters #\n' );
    for n=1:numel(params)
        if params(n).active == 1
            active = 'OPT';
        else
            active = '---';
        end
        value = params(n).value / params(n).step;
        range = params(n).range / params(n).step;
        fprintf( fid, 'description:#%s#:%i:%i:%i:LIN_INT:%s\n', params(n).name, value, range(1), range(2), active );
    end
    fprintf( fid, '#\n\n' );

    fprintf( fid, '# Measurements #\n' );
    fprintf( fid, 'value:---:MIN:0\n' );
    fprintf( fid, '#\n\n' );
    fclose(fid);

    % create extract file
    [a,a,a]=mkdir( 'extract' );
    fid = fopen( 'extract/value', 'wt' );
    fprintf( fid, '# Info #\n' );
    fprintf( fid, 'Name:value\n' );
    fprintf( fid, 'Symbol:value\n' );
    fprintf( fid, '#\n\n' );
    fprintf( fid, '# Commands #\n' );
    fprintf( fid, '#\n\n' );
    fprintf( fid, '# Post Processing #\n' );
    fprintf( fid, 'MEASURE_VAR:   #SYMBOL#: SEARCH_FOR:''value=''\n' );
    fprintf( fid, '#\n' );
    fclose(fid);
    
    % create matlab parameter file
    fid = fopen( 'asco.txt', 'wt' );
    fprintf( fid, '%% this file is processed by asco and variables enclosed in ## are substituted\n' );
    fprintf( fid, 'params = [];\n' );
    for n=1:numel(params)
        fprintf( fid, 'params.%s = #%s# * %i;\n', params(n).name, params(n).name, params(n).step );
    end
    fclose(fid);
    
    % create shell script
    folder_asco_helper = fileparts( mfilename('fullpath') );
    asco_sim_helper = 'optimizer_asco_sim';
    fid = fopen( 'general.sh', 'wt' );
    fprintf( fid, '#!/bin/sh\n' );
    fprintf( fid, 'rm "$2.out" 2> /dev/null\n' );
    fprintf( fid, 'mv "$1.txt" "$1.m"\n' );
    fprintf( fid, 'if [ -f STOP_OPTIMIZATION ]; then\n' );
    fprintf( fid, '    exit\n' );
    fprintf( fid, 'fi\n' );
    fprintf( fid, 'oldpwd=$PWD\n' );
    if isfield(options,'folder_matlabstart')
        % this allows to start the new matlab process in a specific folder
        % => startup.m is picked up here
        fprintf( fid, 'cd "%s"\n', functions.folder_matlabstart );
    end
    if ~isfield(options,'octave_exe')
        % matlab
        fprintf( fid, '%s/bin/matlab -nodesktop -nosplash -r "cd ''%s''; %s(''%s'',''$1'',''$2.out'',''%s''); exit"\n', matlabroot, folder_asco_helper, asco_sim_helper, optimdir, options.simfun );
    else
        % octave
        fprintf( fid, 'export LD_LIBRARY_PATH=\n' );
        fprintf( fid, '%s --silent --eval "cd ''%s''; %s(''%s'',''$1'',''$2.out'',''%s'');"\n', options.octave_exe, folder_asco_helper, asco_sim_helper, optimdir, options.simfun );
    end
    fprintf( fid, 'cd "$oldpwd"\n' );
    fclose(fid);
    fileattrib( 'general.sh', '+x' ); % make it executable

    % clean up old data
    if exist( './best_result.mat', 'file' ), delete( 'best_result.mat' ); end
    if exist( './STOP_OPTIMIZATION', 'file' ), delete( 'STOP_OPTIMIZATION' ); end
    
    % start asco
    [status,result] = unix( 'asco -general asco.txt', '-echo' );
    
    % get best result
    best = load( 'best_result.mat' );
    best = best.best;
    result = best.result;
    for n=1:numel(params)
        name = params(n).name;
        if isfield(best.params,name)
            params(n).value = best.params.(name);
        end
    end
    
elseif strcmp( algorithm, 'simplex-downhill' )
    % ---------------------------------------------------------------------
    % selected algorithm: simplex-downhill
    % Thorsten Liebig <thorsten.liebig@uni-due.de>

    error( 'not implemented yet' );
end

cd( oldfolder );





function folder = absolutepath( folder )
%folder = absolutepath( folder )
% make the path absolute
if isunix
    % Unix
    if folder(1) == '/'
        return
    end
    folder = fullfile( pwd, folder );
else
    % Windows
    folder = strrep( folder, '\', '/' );
    if strcmp( folder(2:3), ':/' ) || strcmp( folder(1:2), '//' ) || (folder(1) == '/')
        return
    end
    if (folder(2) == ':') && (folder(3) ~= '/')
        error( 'relative paths with drive specifier are not supported' );
    end
    folder = fullfile( pwd, folder );
end
