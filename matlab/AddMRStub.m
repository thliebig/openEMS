function CSX = AddMRStub( CSX, materialname, prio, MSL_width, len, alpha, resolution, orientation, normVector, position )
% CSX = AddMRStub( CSX, materialname, prio, MSL_width, len, alpha,
% resolution, orientation, normVector, position )
%
% Microstrip Radial Stub
%
% CSX: CSX-object created by InitCSX()
% materialname: property for the MSL (created by AddMetal() or AddMaterial())
% prio: priority
% MSL_width: width of the MSL to connect the stub to
% len: length of the radial stub
% alpha: angle subtended by the radial stub (degrees) 
% resolution: discrete angle spacing (degrees)
% orientation: angle of main direction of the radial stub (degrees)
% normVector: normal vector of the stub
% position: position of the end of the MSL
%
% This radial stub definition is equivalent to the one Agilent ADS uses.
%
% example:
% CSX = AddMRStub( CSX, 'PEC', 10, 1000, 5900, 30, 1, -90, [0 0 1], [0 -10000 254] );
%
%
% Sebastian Held <sebastian.held@gmx.de>
% Jun 1 2010
%
% See also InitCSX AddMetal AddMaterial

% check normVector
if ~(normVector(1) == normVector(2) == 0) && ...
        ~(normVector(1) == normVector(3) == 0) && ...
        ~(normVector(2) == normVector(3) == 0) || (sum(normVector) == 0)
	error 'normVector must have exactly one component ~= 0'
end
normVector = normVector ./ sum(normVector); % normVector is now a unit vector

% convert angles to radians
alpha_rad = alpha/180*pi;
orientation_rad = orientation/180*pi;
resolution_rad = resolution/180*pi;

%
% build stub at origin (0,0,0) and translate/rotate it later
%

D = 0.5 * MSL_width / sin(alpha_rad/2);
R = cos(alpha_rad/2) * D;

% point at the center of the MSL
p(1,1) = 0;
p(2,1) = -MSL_width/2;
p(1,2) = 0;
p(2,2) = MSL_width/2;

for a = alpha_rad/2 : -resolution_rad : -alpha_rad/2
    p(1,end+1) = cos(a) * (D+len) - R;
    p(2,end)   = sin(a) * (D+len);
end

% rotate
rot = [cos(-orientation_rad), -sin(-orientation_rad); sin(-orientation_rad), cos(-orientation_rad)];
p = (p.' * rot).';

% translate
idx_elevation = [1 2 3];
idx_elevation = idx_elevation(normVector>0);
dim1 = mod( idx_elevation, 3 ) + 1;
dim2 = mod( idx_elevation+1, 3 ) + 1;
p(1,:) = p(1,:) + position(dim1);
p(2,:) = p(2,:) + position(dim2);

elevation = position(idx_elevation);
CSX = AddPolygon( CSX, materialname, prio, idx_elevation-1, elevation, p );
