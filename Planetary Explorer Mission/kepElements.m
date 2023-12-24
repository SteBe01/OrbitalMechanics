%% data - group 2346

addpath("Functions\")
addpath("Functions_custom\")

clear, clc
close all

% orbit data
orbit.a = 0.8016 * 1e4;
orbit.e = 0.1678;
orbit.i = deg2rad(50.3442);
orbit.OM = deg2rad(27.2290);        %Taken from similar object
orbit.om = deg2rad(315.4032);       %Taken from similar object
orbit.theta=deg2rad(122.0796);      %Taken from similar object
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


%% perturbations - cartesian coordinates

n_orbits = 200;
n_points = 10000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan = linspace( 0, T*n_orbits, n_points );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

figure()
plot3( Y(:,1), Y(:,2), Y(:,3))
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag');
axis equal, grid on, hold on
earthPlotSimulation;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or' )
plot3( Y(end,1), Y(end,2), Y(end,3), 'or' )

a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
i_vect = zeros(length(Y), 1);
Om_vect = zeros(length(Y), 1);
om_vect = zeros(length(Y), 1);
theta_vect = zeros(length(Y), 1);

for i = 1:length(Y)
    [a_vect(i), e_vect(i), i_vect(i), Om_vect(i), om_vect(i), theta_vect(i)] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
end

figure()
plot(tspan./(60*60*24), a_vect)
grid on
title('Semi major axis');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24), e_vect)
grid on
title('Eccentricity');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24), rad2deg(i_vect))
grid on
title('Inclination');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24), rad2deg(Om_vect))
grid on
title('RAAN');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(om_vect))
grid on
title('Argument of Periapsis');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(theta_vect))
grid on
title('True Anomaly');
xlabel('time [days]'); ylabel('\theta [°]');


%% perturbations - Gauss's planetary equations

s0 = [orbit.a; orbit.e; orbit.i; orbit.OM; orbit.om; orbit.theta];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, kep] = ode113(@(t,s) eq_motion(t, s, @(t,s) acc_pert_fun(t, s, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, s0, options);

figure()
plot(tspan./(60*60*24), kep(:,1))
grid on
title('Semi major axis');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24), kep(:,2))
grid on
title('Eccentricity');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep(:,3))))
grid on
title('Inclination');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep(:,4))))
grid on
title('RAAN');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep(:,5))))
grid on
title('Argument of Periapsis');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep(:,6))))
grid on
title('True Anomaly');
xlabel('time [days]'); ylabel('\theta [°]');


%% test for movmean

test_a = movmean(kep(:,1), 500);
figure()
plot(test_a)
grid on
title('Semi major axis Movmean Test');
xlabel('time [s]'); ylabel('a [km]');

test_e= movmean(kep(:,2), 500);
figure()
plot(test_e)
grid on
title('Eccentricity Movmean Test');
xlabel('time [s]'); ylabel('e [-]');

test_i = movmean(rad2deg(wrapTo2Pi(kep(:,3))), 500);
figure()
plot(test_i)
grid on
title('Inclination Movmean Test');
xlabel('time [s]'); ylabel('i [°]');

test_Om = movmean(rad2deg(wrapTo2Pi(kep(:,4))), 500);
figure()
plot(test_Om)
grid on
title('RAAN Movmean Test');
xlabel('time [s]'); ylabel('\Omega [°]');

test_om = movmean(rad2deg(wrapTo2Pi(kep(:,5))), 500);
figure()
plot(test_om)
grid on
title('Argument of Periapsis Movmean Test');
xlabel('time [s]'); ylabel('\omega [°]');

test_theta = movmean(rad2deg(wrapTo2Pi(kep(:,6))), 500);
figure()
plot(test_theta)
grid on
title('True Anomaly Movmean Test');
xlabel('time [s]'); ylabel('\theta [°]');


%% All plots together

figure()
plot(tspan./(60*60*24),a_vect)
grid on
hold on 
plot(tspan./(60*60*24),kep(:,1))
plot(tspan./(60*60*24),test_a)
hold off
title('Semi major axis');
legend('Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24),e_vect)
grid on
hold on 
plot(tspan./(60*60*24),kep(:,2))
plot(tspan./(60*60*24),test_e)
hold off
title('Eccentricity');
legend('Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24),rad2deg(i_vect))
grid on
hold on 
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep(:,3))))
plot(tspan./(60*60*24),test_i)
hold off
title('Inclination');
legend('Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24),rad2deg(Om_vect))
grid on
hold on
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep(:,4))))
plot(tspan./(60*60*24),test_Om)
hold off
title('RAAN');
legend('Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24),rad2deg(om_vect))
grid on
hold on 
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep(:,5))))
plot(tspan./(60*60*24),test_om)
hold off
title('Argument of Periapsis');
legend('Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24),rad2deg(theta_vect))
grid on
hold on 
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep(:,6))))
plot(tspan./(60*60*24),test_theta)
hold off
title('True Anomaly');
legend('Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('\theta [°]');





%% other celestial body

addpath("Functions\")
addpath("Functions_custom\")

clear, clc
close all
A = importdata("EXPRESS-MD2.csv");

% orbit data
orbit_new_object.a = A.data(1,10);
orbit_new_object.e = A.data(1,1);
orbit_new_object.i = deg2rad(A.data(1,3));
orbit_new_object.OM = deg2rad(A.data(1,4)); 
orbit_new_object.om = deg2rad(A.data(1,5)); 
orbit_new_object.theta = deg2rad(A.data(1,9)); 
orbit_new_object.kep = [orbit_new_object.a orbit_new_object.e orbit_new_object.i orbit_new_object.OM orbit_new_object.om orbit_new_object.theta];

% Earth data
earth.r = astroConstants(23);
earth.mu = astroConstants(13);
earth.om = deg2rad(15.04) / 3600;
earth.J2 = astroConstants(9);

% perturbation: J2 and Drag (cD = 2.1, A/M = 0.0171 m^2/kg)
spacecraft.cD = 2.1;
spacecraft.AM = 0.0171;


%% without propagation

orbit_new_object.a_no_prop = A.data(:,10);
orbit_new_object.e_no_prop = A.data(:,1);
orbit_new_object.i_no_prop = deg2rad(A.data(:,3));
orbit_new_object.OM_no_prop = deg2rad(A.data(:,4));
orbit_new_object.om_no_prop = deg2rad(A.data(:,5));
orbit_new_object.theta_no_prop = deg2rad(A.data(:,9));

n_orbits = 30;
n_points = length(A.data);

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan= linspace( 0, T*n_orbits, n_points );

figure()
plot(tspan./(60*60*24), orbit_new_object.a_no_prop)
grid on
title('Semi major axis');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24), orbit_new_object.e_no_prop)
grid on
title('Eccentricity');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24), orbit_new_object.i_no_prop)
grid on
title('Inclination');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24), orbit_new_object.OM_no_prop)
grid on
title('RAAN');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24), orbit_new_object.om_no_prop)
grid on
title('Argument of periapsis');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24), orbit_new_object.theta_no_prop)
grid on
title('True Anomaly');
xlabel('time [days]'); ylabel('\theta [°]');


%% perturbations - cartesian coordinates

kep_body = [A.data(1,10), A.data(1,1), deg2rad(A.data(1,3)), deg2rad(A.data(1,4)), deg2rad(A.data(1,5)), deg2rad(A.data(1,9)), earth.mu];
[r0, v0] = kep2car(kep_body);

n_orbits = 30;
n_points = length(A.data);

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan= linspace( 0, T*n_orbits, n_points );

y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

% Plot orbit perturbed
figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag');
axis equal, grid on, hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or' )
plot3( Y(end,1), Y(end,2), Y(end,3), 'or' )

a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
i_vect = zeros(length(Y), 1);
Om_vect = zeros(length(Y), 1);
om_vect = zeros(length(Y), 1);
theta_vect = zeros(length(Y), 1);

for i = 1:length(Y)
    [a_vect(i), e_vect(i), i_vect(i), Om_vect(i), om_vect(i), theta_vect(i)] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
end

figure()
plot(tspan./(60*60*24), a_vect)
grid on
title('Semi major axis');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24), e_vect)
grid on
title('Eccentricity');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24), rad2deg(i_vect))
grid on
title('Inclination');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24), rad2deg(Om_vect))
grid on
title('RAAN');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(om_vect))
grid on
title('Argument of periapsis');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(theta_vect))
grid on
title('True Anomaly');
xlabel('time [days]'); ylabel('\theta [°]');


%% perturbations - Gauss's planetary equations

s0 = [orbit_new_object.a; orbit_new_object.e; orbit_new_object.i; orbit_new_object.OM; orbit_new_object.om; orbit_new_object.theta];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, kep_new_object] = ode113(@(t,s) eq_motion(t, s, @(t,s) acc_pert_fun(t, s, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, s0, options);

figure()
plot(tspan./(60*60*24), kep_new_object(:,1))
grid on
title('Semi major axis');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24), kep_new_object(:,2))
grid on
title('Eccentricity');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep_new_object(:,3))))
grid on
title('Inclination');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep_new_object(:,4))))
grid on
title('RAAN');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep_new_object(:,5))))
grid on
title('Argument of Periapsis');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24), rad2deg(wrapTo2Pi(kep_new_object(:,6))))
grid on
title('True Anomaly');
xlabel('time [days]'); ylabel('\theta [°]');


%% test for movmean

test_a = movmean(kep_new_object(:,1), 500);
figure()
plot(tspan./(60*60*24),test_a)
grid on
title('Semi major axis Movmean Test');
xlabel('time [days]'); ylabel('a [km]');

test_e= movmean(kep_new_object(:,2), 500);
figure()
plot(tspan./(60*60*24),test_e)
grid on
title('Eccentricity Movmean Test');
xlabel('time [days]'); ylabel('e [-]');

test_i = movmean(rad2deg(wrapTo2Pi(kep_new_object(:,3))), 500);
figure()
plot(tspan./(60*60*24),test_i)
grid on
title('Inclination Movmean Test');
xlabel('time [days]'); ylabel('i [°]');

test_Om = movmean(rad2deg(wrapTo2Pi(kep_new_object(:,4))), 500);
figure()
plot(tspan./(60*60*24),test_Om)
grid on
title('RAAN Movmean Test');
xlabel('time [days]'); ylabel('\Omega [°]');

test_om = movmean(rad2deg(wrapTo2Pi(kep_new_object(:,5))), 500);
figure()
plot(tspan./(60*60*24),test_om)
grid on
title('Argument of Periapsis Movmean Test');
xlabel('time [days]'); ylabel('\omega [°]');

test_theta = movmean(rad2deg(wrapTo2Pi(kep_new_object(:,6))), 500);
figure()
plot(tspan./(60*60*24),test_theta)
grid on
title('True Anomaly Movmean Test');
xlabel('time [days]'); ylabel('\theta [°]');


%% All plots together

figure()
plot(tspan./(60*60*24), orbit_new_object.a_no_prop)
grid on
hold on
plot(tspan./(60*60*24),a_vect)
plot(tspan./(60*60*24),kep_new_object(:,1))
plot(tspan./(60*60*24),test_a)
hold off
title('Semi major axis');
legend('Real Data','Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('a [km]');

figure()
plot(tspan./(60*60*24), orbit_new_object.e_no_prop)
grid on
hold on 
plot(tspan./(60*60*24),e_vect)
plot(tspan./(60*60*24),kep_new_object(:,2))
plot(tspan./(60*60*24),test_e)
hold off
title('Eccentricity');
legend('Real Data','Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('e [-]');

figure()
plot(tspan./(60*60*24), orbit_new_object.i_no_prop)
grid on
hold on
plot(tspan./(60*60*24),rad2deg(i_vect))
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep_new_object(:,3))))
plot(tspan./(60*60*24),test_i)
hold off
title('Inclination');
legend('Real Data','Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('i [°]');

figure()
plot(tspan./(60*60*24), orbit_new_object.OM_no_prop)
grid on
hold on
plot(tspan./(60*60*24),rad2deg(Om_vect))
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep_new_object(:,4))))
plot(tspan./(60*60*24),test_Om)
hold off
title('RAAN');
legend('Real Data','Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('\Omega [°]');

figure()
plot(tspan./(60*60*24), orbit_new_object.om_no_prop)
grid on
hold on 
plot(tspan./(60*60*24),rad2deg(om_vect))
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep_new_object(:,5))))
plot(tspan./(60*60*24),test_om)
hold off
title('Argument of Periapsis');
legend('Real Data','Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('\omega [°]');

figure()
plot(tspan./(60*60*24), orbit_new_object.theta_no_prop)
grid on
hold on 
plot(tspan./(60*60*24),rad2deg(theta_vect))
plot(tspan./(60*60*24),rad2deg(wrapTo2Pi(kep_new_object(:,6))))
plot(tspan./(60*60*24),test_theta)
hold off
title('True Anomaly');
legend('Real Data','Cartesian','Gauss','Filtered');
xlabel('time [days]'); ylabel('\theta [°]');


%% movmean finder

% data = A.data(:,1);
% 
% figure
% hold on
% 
% for i = 1:1:100
%     a = plot(movmean(data, i));
%     drawnow
%     i
%     pause(1)
%     delete(a)
% end
% legend

% 33, 65, 98
% 19, 25, 31, 37, 50

