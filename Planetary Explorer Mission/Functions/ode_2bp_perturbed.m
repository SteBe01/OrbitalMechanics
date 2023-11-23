function dy = ode_2bp_perturbed( ~, y, mu, Re, J2, om_E, A_M, cD)
%ode_2bp ODE system for the two-body problem (Keplerian motion)
%
% PROTOTYPE
% dy = ode_2bp( t, y, mu, Re, J2 )
%
% INPUT:
% t[1] Time (can be omitted, as the system is autonomous) [T]
% y[6x1] State of the body ( rx, ry, rz, vx, vy, vz ) [ L, L/T ]
% mu[1] Gravitational parameter of the primary [L^3/T^2]
% Re[1] Earth’s equatorial radius
% J2[1] constant coefficient
%
% OUTPUT:
% dy[6x1] Derivative of the state [ L/T^2, L/T^3 ]
%
% CONTRIBUTORS:
% Juan Luis Gonzalo Gomez
%
% VERSIONS
% 2018-09-26: First version
%
% -------------------------------------------------------------------------

% Position and velocity
r = y(1:3);
v = y(4:6);

% Distance from the primary
rnorm = norm(r);

% Calculate perturbation
perturbCoeff = 3/2 * (J2*mu*Re^2)/rnorm^4;
a_j2 = perturbCoeff * [r(1)/rnorm*(5*r(3)^2/rnorm^2-1); ...
                       r(2)/rnorm*(5*r(3)^2/rnorm^2-1); ...
                       r(3)/rnorm*(5*r(3)^2/rnorm^2-3)];

rho0 = 3.019 * 1e-15;
h0 = 1000;
H = 268;
h = norm(r) - Re;
rho = rho0 * exp(-(h-h0)/H);
v_rel = v - cross([0 0 om_E], r)';

v_rel = v_rel*1000;
a_drag = -0.5 * A_M * rho * cD * norm(v_rel)^2 * (v_rel ./ norm(v_rel));

% Set the derivatives of the state
dy = [ v
 (-mu/rnorm^3)*r + a_j2 + a_drag];
end