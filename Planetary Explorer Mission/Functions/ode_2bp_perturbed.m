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

% J2
perturbCoeff = 3/2 * (J2*mu*Re^2)/rnorm^4;
a_j2 = perturbCoeff * [r(1)/rnorm*(5*r(3)^2/rnorm^2-1); ...
                       r(2)/rnorm*(5*r(3)^2/rnorm^2-1); ...
                       r(3)/rnorm*(5*r(3)^2/rnorm^2-3)];


% air drag
h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 ...
    300 350 400 450 500 600 700 800 900 1000]';
rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3, 3.206*1e-4 ...
    8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 2.438*1e-8 ...
    8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 2.789*1e-10 7.248*1e-11 ...
    2.418*1e-11 9.158*1e-12 3.725*1e-12 1.585*1e-12 6.967*1e-13 1.454*1e-13...
    3.614*1e-14 1.170*1e-14 5.245*1e-15 3.019*1e-15]';
H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
    7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 53.298 ...
    58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

h = norm(r) - Re;

found = 0;
for i = 1:length(h0_vect)-1
    if h0_vect(i+1) > h
        found = 1;
        break
    end
end
if ~found
    i = length(h0_vect);
end

rho0 = rho0_vect(i);
h0 = h0_vect(i);
H = H_vect(i);

rho = rho0 * exp(-(h-h0)/H);
v_rel = v - cross([0 0 om_E], r)';
v_rel = v_rel * 1e3;

a_drag = (-0.5 * A_M * rho * cD * norm(v_rel)^2) .* (v_rel ./ norm(v_rel)) .* 1e-3;

% Set the derivatives of the state
dy = [ v
 (-mu/rnorm^3)*r + a_j2 + a_drag];

if h <= 8.848
    error("Impact on Earth detected! Reduce simulation time");
end


end
