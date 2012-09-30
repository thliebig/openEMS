function [E_theta,E_phi,Prad,Dmax] = AnalyzeNF2FF( Sim_Path, nf2ff, f, theta, phi, r )
% [E_theta,E_phi,Prad,Dmax] = AnalyzeNF2FF( Sim_Path, filenames_E, filenames_H, f, theta, phi, r )
%
% calculates the farfield via a near field to far field transformation
%
% input:
%   Sim_Path:    simulation directory
%   nf2ff:       structure on filenames etc. as created by CreateNF2FFBox 
%   f:           frequency (Hz) for far field calculation
%   theta:       (degrees) vector of discrete theta values to calculate the far field for
%   phi:         (degrees) vector of discrete phi values to calculate the far field for
%   r:           (optional) Radius (m) at which the E-fields are calculated (default: 1 m)
%
% output:
%   E_theta:     E_theta(theta,phi); theta component of the electric field strength at radius r
%   E_phi:       E_phi(theta,phi);     phi component of the electric field strength at radius r
%   Prad:        time averaged radiated power
%   Dmax:        maximum directivity
%
% example:
%   see examples/NF2FF/infDipol.m
%
% See also CreateNF2FFBox
%
% (C) 2010 Sebastian Held <sebastian.held@gmx.de>
% (C) 2011 Thorsten Liebig <thorsten.liebig@gmx.de>

% check arguments
error( nargchk(6,6,nargin) );
if ~isscalar(f)
    error 'Currently only one frequency is supported. Call this function multiple times.'
end

warning('openEMS:AnalyzeNF2FF','This function is deprecated, use CalcNF2FF instead');

filenames_E = nf2ff.filenames_E;
filenames_H = nf2ff.filenames_H;

if (~isrow(nf2ff.directions))
    nf2ff.directions = nf2ff.directions';
end

% read time domain field data and transform into frequency domain
for n=find(nf2ff.directions==1)
    [Ef{n}, E_mesh{n}] = ReadHDF5Dump( [Sim_Path '/' filenames_E{n} '.h5'], 'Frequency', f );

    if (Ef{n}.FD.frequency(1) ~= f)
        error 'frequency mismach'
    end

    %clear out time domain data
    if isfield(Ef{n},'TD')
        Ef{n} = rmfield(Ef{n},'TD');
    end

    [Hf{n}, H_mesh{n}] = ReadHDF5Dump( [Sim_Path '/' filenames_H{n} '.h5'], 'Frequency', f );
    %clear out time domain data
    if isfield(Hf{n},'TD')
        Hf{n} = rmfield(Hf{n},'TD');
    end

    % reshape mesh into row vector
    mesh{n}.x = reshape( E_mesh{n}.lines{1}, 1, [] );
    mesh{n}.y = reshape( E_mesh{n}.lines{2}, 1, [] );
    mesh{n}.z = reshape( E_mesh{n}.lines{3}, 1, [] );
end

% create a normal vector for every plane
% FIXME!!! this is dependent upon the order of filenames_*
n = {};
for a=1:6
    temp = [(a<=2), ((a>=3)&&(a<=4)), (a>=5)];
    n{a} = temp - 2*mod(a,2)*temp;
end


physical_constants

k = 2*pi*f/c0;
center = [0 0 0];
Umax = 0;

phi_idx = 0;
for phi_deg_aufpunkt = phi
    phi_rad_aufpunkt = phi_deg_aufpunkt/180*pi; % radiant
    phi_idx = phi_idx + 1;
    
    theta_idx = 0;
    for theta_deg_aufpunkt = theta
        theta_rad_aufpunkt = theta_deg_aufpunkt/180*pi; % radiant
        theta_idx = theta_idx + 1;

        N_theta = 0;
        N_phi = 0;
        L_theta = 0;
        L_phi = 0;
        for a=find(nf2ff.directions==1)
            [N_theta_,N_phi_,L_theta_,L_phi_] = process_plane( k, n{a}, center, mesh{a}, Ef{a}.FD.values{1}, Hf{a}.FD.values{1}, theta_rad_aufpunkt, phi_rad_aufpunkt );
            N_theta = N_theta + N_theta_; N_phi = N_phi + N_phi_;
            L_theta = L_theta + L_theta_; L_phi = L_phi + L_phi_;
        end

        % E-fields
        erg_E_theta = -1i*k*exp(-1i*k*r) / (4*pi*r)*(L_phi+Z0*N_theta);
        erg_E_phi   =  1i*k*exp(-1i*k*r) / (4*pi*r)*(L_theta-Z0*N_phi);

        % output
        E_theta(theta_idx,phi_idx) = erg_E_theta;
        E_phi(theta_idx,phi_idx)   = erg_E_phi;

        % directivity
        U = r^2/(2*Z0) * sum(abs([erg_E_theta erg_E_phi]).^2);
        Umax = max( [Umax U] );
    end
end

% power
Prad = 0;
for a=find(nf2ff.directions==1)
    [~,~,~,~,P] = process_plane( k, n{a}, center, mesh{a}, Ef{a}.FD.values{1}, Hf{a}.FD.values{1}, theta_rad_aufpunkt, phi_rad_aufpunkt );
    Prad = Prad + P;
end

% directivity
Dmax = 4*pi*Umax / Prad;


% integrate over one plane
function [N_theta,N_phi,L_theta,L_phi,Prad] = process_plane( k, n, center, mesh, E_field, H_field, theta_rad_aufpunkt, phi_rad_aufpunkt )
% [N_theta,N_phi,L_theta,L_phi,Prad] = process_plane( k, n, center, mesh, E_field, H_field, theta_rad_aufpunkt, phi_rad_aufpunkt )
%
% k: wave number
% n: normal vector of the plane
% center: correction coordinates for the center of the antenna
% mesh: mesh info
% E_field: E field array ?x?x?x3
% H_field: H field array ?x?x?x3

% speed up
sin__theta_rad_aufpunkt = sin(theta_rad_aufpunkt);
cos__theta_rad_aufpunkt = cos(theta_rad_aufpunkt);
sin__phi_rad_aufpunkt   = sin(phi_rad_aufpunkt);
cos__phi_rad_aufpunkt   = cos(phi_rad_aufpunkt);

if abs(n(1)) == 1
    % x-plane
    x = mesh.x(1);
    [y z] = ndgrid( mesh.y, mesh.z );
    coord1 = mesh.y.';
    coord2 = mesh.z.';
    Ex = squeeze( E_field(1,:,:,1) );
    Ey = squeeze( E_field(1,:,:,2) );
    Ez = squeeze( E_field(1,:,:,3) );
    Hx = squeeze( H_field(1,:,:,1) );
    Hy = squeeze( H_field(1,:,:,2) );
    Hz = squeeze( H_field(1,:,:,3) );
elseif abs(n(2)) == 1
    % y-plane
    y = mesh.y(1);
    [x z] = ndgrid( mesh.x, mesh.z );
    coord1 = mesh.x.';
    coord2 = mesh.z.';
    Ex = squeeze( E_field(:,1,:,1) );
    Ey = squeeze( E_field(:,1,:,2) );
    Ez = squeeze( E_field(:,1,:,3) );
    Hx = squeeze( H_field(:,1,:,1) );
    Hy = squeeze( H_field(:,1,:,2) );
    Hz = squeeze( H_field(:,1,:,3) );
elseif abs(n(3)) == 1
    % z-plane
    z = mesh.z(1);
    [x y] = ndgrid( mesh.x, mesh.y );
    coord1 = mesh.x.';
    coord2 = mesh.y.';
    Ex = squeeze( E_field(:,:,1,1) );
    Ey = squeeze( E_field(:,:,1,2) );
    Ez = squeeze( E_field(:,:,1,3) );
    Hx = squeeze( H_field(:,:,1,1) );
    Hy = squeeze( H_field(:,:,1,2) );
    Hz = squeeze( H_field(:,:,1,3) );
end    

Jx = n(2) .* Hz - n(3) .* Hy;
Jy = n(3) .* Hx - n(1) .* Hz;
Jz = n(1) .* Hy - n(2) .* Hx;
Mx = -n(2) .* Ez + n(3) .* Ey;
My = -n(3) .* Ex + n(1) .* Ez;
Mz = -n(1) .* Ey + n(2) .* Ex;
r_cos_psi = x*sin__theta_rad_aufpunkt*cos__phi_rad_aufpunkt + y*sin__theta_rad_aufpunkt*sin__phi_rad_aufpunkt + z*cos__theta_rad_aufpunkt;
e_fkt     = exp( +1i*k*r_cos_psi );
N_theta = dbltrapz( ( Jx*cos__theta_rad_aufpunkt*cos__phi_rad_aufpunkt + Jy*cos__theta_rad_aufpunkt*sin__phi_rad_aufpunkt - Jz*sin__theta_rad_aufpunkt) .* e_fkt, coord1, coord2 );
N_phi   = dbltrapz( (-Jx*sin__phi_rad_aufpunkt + Jy*cos__phi_rad_aufpunkt) .* e_fkt, coord1, coord2 );
L_theta = dbltrapz( ( Mx*cos__theta_rad_aufpunkt*cos__phi_rad_aufpunkt + My*cos__theta_rad_aufpunkt*sin__phi_rad_aufpunkt - Mz*sin__theta_rad_aufpunkt) .* e_fkt, coord1, coord2 );
L_phi   = dbltrapz( (-Mx*sin__phi_rad_aufpunkt + My*cos__phi_rad_aufpunkt) .* e_fkt, coord1, coord2 );

if nargout > 4
    % Prad requested
    
    % this is crap! recode it!
    EH = zeros(size(Ex));
    for i1 = 1:numel(coord1)
        for i2 = 1:numel(coord2)
            E = [Ex(i1,i2) Ey(i1,i2) Ez(i1,i2)];
            H = [Hx(i1,i2) Hy(i1,i2) Hz(i1,i2)];
            EH(i1,i2) = real( dot(cross(E,conj(H)),n) );
        end
    end
    Prad = 0.5 * dbltrapz( EH, coord1, coord2 );
end




function Q = dbltrapz(matrix,a,b)
%DBLTRAPZ  Trapezoidal numerical integration in two dimensions.
%   Z = DBLTRAPZ(MATRIX,A,B) computes an approximation of the double integral
%   of MATRIX via the trapezoidal method (with respect to A and B).  A and B must be
%   column vectors of the same length.
%   index like this: MATRIX(A,B)

if nargin < 3, error('MATLAB:dblquad:NotEnoughInputs',...
                     'Requires at least three inputs.'); end
if size(a,2) ~= 1, error('column vectors required'); end
if size(b,2) ~= 1, error('column vectors required'); end

temp = zeros(size(b));
for i = 1:length(b) 
    temp(i) = trapz( a, matrix(:,i) );
end

Q = trapz( b, temp );
