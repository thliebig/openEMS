function result = optimizer_simfun( folder, params )
%
% simulation function
%
% the variable params contains the simulation parameters

disp( [mfilename ': SIMULATING...'] );

if nargin == 0
    % visualize the structure if called without parameters
    folder = 'tmp';
    params.turns = 10;
end

oldpwd = pwd;
[a,a,a] = mkdir( folder );
cd( folder );

% create the structure
f_max = 50e6;
structure( params, 'tmp.xml', f_max );

if nargin == 0
    % visualize the structure
    CSXGeomPlot('tmp.xml');
    return;
end

% start simulation
RunOpenEMS( '.', 'tmp.xml', '--engine=fastest' );

% postprocess the results
L = postproc( 'tmp.xml', f_max );
disp( ['DONE. L = ' num2str(L(1)/1e-6) ' uH'] );

% calculate result
goal   = 2e-6; % specify the goal: 2 uH
result = abs(goal - L(1)); % costs must not be negative

% restore curent folder
cd( oldpwd );




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

function CSX = structure( params, filename, f_max )
% CSX = structure( params, filename )

unit   = 1e-3; % specify length in mm
lambda = 3e8/f_max;
resolution = lambda/15/unit;
mesh_size  = 1.5;
radius = 10;
turns  = params.turns;
height = 4 * turns;
feed_length = 10;

CSX = InitCSX();

%% create coil
p1 = create_coil( radius, height, turns );
p  = p1(:,end) + [feed_length;0;0];
p1 = [p1 p];
p  = p1(:,1) + [feed_length;0;0];
p1 = [p p1];
CSX = AddMetal(CSX,'PEC1');
CSX = AddCurve(CSX, 'PEC1', 0, p1);

%% create mesh
extraspace = 5*radius;
mesh.x = linspace(-radius,radius,ceil(2*radius/mesh_size));
mesh.x = [mesh.x mesh.x(1)-extraspace mesh.x(end)+extraspace];
mesh.x = [mesh.x p1(1,1) p1(1,1)-mesh_size p1(1,1)+mesh_size];
mesh.x = SmoothMeshLines2( mesh.x, resolution );
mesh.y = linspace(-radius,radius,ceil(2*radius/mesh_size));
mesh.y = [mesh.y mesh.y(1)-extraspace mesh.y(end)+extraspace];
mesh.y = SmoothMeshLines2( mesh.y, resolution );
mesh.z = linspace(0,height,ceil(height/mesh_size));
mesh.z = [mesh.z mesh.z(1)-extraspace mesh.z(end)+extraspace];
% mesh.z = [mesh.z p1(3,1) p1(3,1)-mesh_size p1(3,1)+mesh_size];
mesh.z = SmoothMeshLines2( mesh.z, resolution );
CSX = DefineRectGrid(CSX, unit, mesh);


%% create port
[CSX,port] = AddCurvePort( CSX, 10, 1, 50, p1(:,1), p1(:,end), 'excite' );

if nargin > 1
    max_timesteps = 100000;
    min_decrement = 1e-5;
    FDTD = InitFDTD( max_timesteps, min_decrement );
    FDTD = SetGaussExcite( FDTD, 0, f_max );
    BC = {'PEC' 'PEC' 'PEC' 'PEC' 'PMC' 'PMC'};
    FDTD = SetBoundaryCond( FDTD, BC );

    WriteOpenEMS( filename, FDTD, CSX );
end


function p = create_coil(coil_rad,coil_length,coil_turns,coil_res,winding_direction,direction,offset,angle_offset)
if nargin < 8, angle_offset = 0; end
if nargin < 7, offset = [0; 0; 0]; end
if nargin < 6, direction = +1; end
if nargin < 5, winding_direction = +1; end
if nargin < 4, coil_res = 30; end
dt = 1/coil_res;
height = 0;

p = [];
while abs(height) < coil_length
    angle      = height / (coil_length/coil_turns) * 2*pi;
    p(1,end+1) = coil_rad * cos(angle*winding_direction+angle_offset);
    p(2,end)   = coil_rad * sin(angle*winding_direction+angle_offset);
    p(3,end)   = height * direction;
    p(:,end)   = p(:,end) + offset;
    height     = height + coil_length/coil_turns * dt;
end



function L = postproc( filename, f_max )
freq = linspace(0,f_max,201);
freq(1) = []; % delete DC component

folder = fileparts(filename);
U = ReadUI( 'port_ut1', folder, freq );
I = ReadUI( 'port_it1', folder, freq );
Z = U.FD{1}.val ./ I.FD{1}.val;

L = imag(Z) ./ (2*pi*freq);
L = reshape( L, 1, [] ); % row vector
