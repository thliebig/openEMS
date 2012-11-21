%
% run the testsuite
%

clc
clear
close all
drawnow

if isOctave
    confirm_recursive_rmdir(0);
    page_screen_output(0);      % do not buffer output
    page_output_immediately(1); % do not buffer output
end

folder = fileparts( mfilename( 'fullpath' ) );
cd( folder );
addpath( [folder filesep 'helperscripts'] );

% openEMS options
options = {'--engine=multithreaded', '--engine=sse-compressed', '--engine=sse', '--engine=basic'};

for o=1:numel(options)

    disp( [datestr(now) ' *** TESTSUITE started (options: ' options{o} ')'] );

    % now list the tests
    folders = dir();
    for f=1:numel(folders)
        if folders(f).isdir
            if strcmp(folders(f).name,'.') || strcmp(folders(f).name,'..')
                continue
            end
            if strcmp(folders(f).name,'helperscripts')
                continue
            end
            oldpwd = pwd;
            cd( folders(f).name );
            scripts = dir('*.m');
            for s=1:numel(scripts)
                if ~scripts(s).isdir
                    % execute function
                    disp( [datestr(now) ' executing: ' folders(f).name '/' scripts(s).name] );
                    [~,fname] = fileparts( scripts(s).name );
                    if isOctave
                        fflush(1); % flush stdout
                    end
                    pass = feval( fname, options{o}, 'run_testsuite' );
                end
            end
            cd(oldpwd);
        end
    end
end

disp( '***' );
disp( ['*** ' datestr(now) ' ALL TESTS DONE'] );
disp( '***' );
