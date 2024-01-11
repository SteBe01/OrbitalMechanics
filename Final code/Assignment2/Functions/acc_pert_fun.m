function acc_pert_vec = acc_pert_fun(t, y, mu, Re, J2, om_E, A_M, cD)

% Calculates perturbing accelerations
%
% Usage:
% acc_pert_vec = acc_pert_fun(t, y, mu, Re, J2, om_E, A_M, cD)
%
% Input arguments:
% ----------------------------------------------------------------
% t             [Nx1]   time                        [-/s]
% y             [1x6]   Keplerian Elements Vector   [-]
% mu            [1x1]   Planetary constant          [km^3/s^2]
% Re            [1x1]   Mean radius of the planet   [km]
% J2            [1x1]   Gravitatonal field constant [-]
% om_E          [1x1]   Earth angular velocity      [rad/s]
% A_M           [1x1]   Reference area over mass    [m^2/kg]
% cD            [1x1]   Drag coefficient            [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% acc_pert_vec  [1x3]   Perturbing Accelerations    [km/s^2]
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

a = y(1);
e = y(2);
i = y(3);
OM = y(4);
om = y(5);
th = y(6);

p = a*(1-e^2);
r = p/(1+e*cos(th));
h = sqrt(p*mu);
v = sqrt(2*mu/r - mu/a);

[rr, vv] = kep2car(a, e, i, OM, om, th, mu);

perturbCoeff = 3/2 * (J2*mu*Re^2)/r^4;
a_j2_xyz = perturbCoeff * [rr(1)/r*(5*rr(3)^2/r^2-1); ...
                   rr(2)/r*(5*rr(3)^2/r^2-1); ...
                   rr(3)/r*(5*rr(3)^2/r^2-3)];

u = th + om;
Qxr = [-sin(OM)*cos(i)*sin(u) + cos(OM)*cos(u)   cos(OM)*cos(i)*sin(u) + sin(OM)*cos(u)  sin(i)*sin(u)
       -sin(OM)*cos(i)*cos(u) - cos(OM)*sin(u)   cos(OM)*cos(i)*cos(u) - sin(OM)*sin(u)  sin(i)*cos(u)
       sin(OM)*sin(i)                           -cos(OM)*sin(i)                         cos(i)];

a_j2_rsw = Qxr*a_j2_xyz;

rot_mat = h/(p*v) * [e*sin(th)      -(1+e*cos(th));
                     1+e*cos(th)    e*sin(th)];

a_j2_tn = rot_mat\a_j2_rsw(1:2);
a_j2_tnh = [a_j2_tn; a_j2_rsw(3)];

% Air drag
h0_vect = [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 ...
    180 200 250 300 350 400 450 500 600 700 800 900 1000]';
rho0_vect = [1.225 3.899*1e-2 1.774*1e-2 3.972*1e-3 1.057*1e-3 ...
    3.206*1e-4 8.770*1e-5 1.905*1e-5 3.396*1e-6 5.297*1e-7 9.661*1e-8 ...
    2.438*1e-8 8.484*1e-9 3.845*1e-9 2.070*1e-9 5.464*1e-10 ...
    2.789*1e-10 7.248*1e-11 2.418*1e-11 9.158*1e-12 3.725*1e-12 ...
    1.585*1e-12 6.967*1e-13 1.454*1e-13 3.614*1e-14 1.170*1e-14 ...
    5.245*1e-15 3.019*1e-15]';
H_vect = [7.249 6.349 6.682 7.554 8.382 7.714 6.549 5.799 5.382 5.877 ...
    7.263 9.473 12.636 16.149 22.523 29.740 37.105 45.546 53.628 ...
    53.298 58.515 60.828 63.822 71.835 88.667 124.64 181.05 268.00]';

height = norm(rr) - Re;                         % height from ground

idx = sum(h0_vect < height);

rho0 = rho0_vect(idx);
h0 = h0_vect(idx);                              % [km]
H = H_vect(idx);                                % [km]

rho = rho0 * exp(-(height-h0)/H);               % [kg/m^3]

v_rel = vv - cross([0 0 om_E], rr');            % [km]
v_rel = v_rel * 1e3;                            % [m]

a_drag_xyz = (-0.5 * A_M * rho * cD * norm(v_rel)^2) ...
    * (v_rel ./ norm(v_rel)) * 1e-3;            % [m/s^2]]

a_drag_rsw = Qxr*a_drag_xyz';

a_drag_tn = rot_mat\a_drag_rsw(1:2);
a_drag_tnh = [a_drag_tn; a_drag_rsw(3)];

if height <= 0
    disp("current_time: " + num2str(t) + " parameters: " + num2str(y));
    error("Impact on ground detected");
end

acc_pert_vec = a_j2_tnh + a_drag_tnh;

end

