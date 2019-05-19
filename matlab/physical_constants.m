%
% physical constants
%

% Bronstein 3rd ed., 1997, pp. 945-946
C0 = 299792458; % m/s
c0 = C0; % constants in capital letters, c0 for legacy support
MUE0 = 4e-7*pi; % N/A^2
EPS0 = 1/(MUE0*C0^2); % F/m

% free space wave impedance
Z0 = sqrt(MUE0/EPS0); % Ohm
