%% data - group 2346

addpath("Functions\")
addpath("Functions_custom\")

clear, clc
close all

% orbit data
orbit.a = 0.8016 * 1e4;
orbit.e = 0.1678;
orbit.i = deg2rad(50.3442);
orbit.OM = deg2rad(27.2290); %Taken from similar object
orbit.om = deg2rad(315.4032); %Taken from similar object
orbit.theta=deg2rad(122.0796); %Taken from similar object
orbit.kep = [orbit.a orbit.e orbit.i orbit.OM orbit.om orbit.theta];
orbit.ratio_k = 12;
orbit.ratio_m = 1;

% Earth data
earth.r = astroConstants(23);
earth.mu = astroConstants(13);
earth.om = deg2rad(15.04) / 3600;
earth.J2 = astroConstants(9);

% perturbation: J2 and Drag (cD = 2.1, A/M = 0.0171 m^2/kg)
spacecraft.cD = 2.1;
spacecraft.AM = 0.0171;


%% J2 Secular Effect

T = 2*pi*sqrt( orbit.a^3/earth.mu )*100;
tspan= linspace( 0, T, 100 );

[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,earth.mu), tspan, y0, options );

a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
i_vect = zeros(length(Y), 1);
Om_vect = zeros(length(Y), 1);
om_vect = zeros(length(Y), 1);
theta_vect = zeros(length(Y), 1);
Om_sec = zeros(length(Y), 1);
om_sec = zeros(length(Y), 1);
theta_sec = zeros(length(Y), 1);

for i = 1:length(Y)
    [a_vect(i), e_vect(i), i_vect(i), Om_vect(i), om_vect(i), theta_vect(i)] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
    parameter = -((3/2)*(sqrt(earth.mu)*earth.J2*earth.r^2)/(((1-e_vect(i)^2)^2))*a_vect(i)^(7/2));
    parameter_2 = ((3/2)*(sqrt(earth.mu)*earth.J2*earth.r^2)/(((1-e_vect(i)^(3/2))^2))*a_vect(i)^(7/2));
    Om_sec(i) = parameter*cos(i_vect(i)); 
    om_sec(i) = parameter*(((5/2)*(sin(i_vect(i)))^2)-2); 
    theta_sec(i) = parameter_2*(1-((3/2)*((sin(i_vect(i)))^2))); 

end

figure()
plot(rad2deg(Om_sec))
grid on
title('RAAN');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(om_sec))
grid on
title('Argument of periapsis');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(theta_sec))
grid on
title('True Anomaly');
xlabel('time [s]'); ylabel('\theta [°]');