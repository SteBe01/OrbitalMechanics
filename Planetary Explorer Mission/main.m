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


%% unperturbed 2bp orbit

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, 100 );

[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,earth.mu), tspan, y0, options );

figure()
comet3( Y(:,1), Y(:,2), Y(:,3))
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

title('Two-body problem orbit');
axis equal, grid on,hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or' )


%% ground track plot (unperturbed)

om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)

orbit_number = orbit.ratio_k*10;
tspan_dim = 100000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot(lon, lat, "EarthTexture.jpg")


orbit.a_rep = aFinder(orbit.ratio_k, orbit.ratio_m, om_E, earth.mu);
T = 2*pi*sqrt( orbit.a_rep^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car(orbit.a_rep, orbit.e, orbit.i, orbit.OM, orbit.om, 0, earth.mu);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot(lon, lat, "EarthTexture.jpg")


%% perturbations

n_orbits = 200;
n_points = 10000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T*n_orbits, n_points );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, 0, earth.mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

figure()
comet3( Y(:,1), Y(:,2), Y(:,3))
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
plot(a_vect)
grid on
title('Semi major axis');
xlabel('time [s]'); ylabel('a [km]');

figure()
plot(e_vect)
grid on
title('Eccentricity');
xlabel('time [s]'); ylabel('e [-]');

figure()
plot(rad2deg(i_vect))
grid on
title('Inclination');
xlabel('time [s]'); ylabel('i [°]');

figure()
plot(rad2deg(Om_vect))
grid on
title('RAAN');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(om_vect))
grid on
title('Argument of Periapsis');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(theta_vect))
grid on
title('True Anomaly');
xlabel('time [s]'); ylabel('\theta [°]');

%% J2 Secular Effect

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
    Om_sec = parameter*cos(i_vect(i)); 
    om_sec =parameter*(((5/2)*(sin(i_vect(i)))^2)-2); 
    theta_sec = parameter_2*(1-((2/3)*((sin(i_vect(i)))^2))); 

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


%% Perturbation over time

kep0 = [orbit.a; orbit.e; orbit.i; orbit.OM; orbit.om; 0];
[t, kep] = ode113(@(t,kep) Pert_guass_eq_tnh_frame(t, kep, @(t,kep) acc_pert_fun(t, kep, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, kep0, options);

figure()
plot(kep(:,1))
grid on
title('Semi major axis');
xlabel('time [s]'); ylabel('a [km]');

figure()
plot(kep(:,2))
grid on
title('Eccentricity');
xlabel('time [s]'); ylabel('e [-]');

figure()
plot(rad2deg(kep(:,3)))
grid on
title('Inclination');
xlabel('time [s]'); ylabel('i [°]');

figure()
plot(rad2deg(kep(:,4)))
grid on
title('RAAN');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(kep(:,5)))
grid on
title('Argument of Periapsis');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(kep(:,6)))
grid on
title('True Anomaly');
xlabel('time [s]'); ylabel('\theta [°]');


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

test_i = movmean(rad2deg(kep(:,3)), 500);
figure()
plot(test_i)
grid on
title('Inclination Movmean Test');
xlabel('time [s]'); ylabel('i [°]');

test_Om = movmean(rad2deg(kep(:,4)), 500);
figure()
plot(test_Om)
grid on
title('RAAN Movmean Test');
xlabel('time [s]'); ylabel('\Omega [°]');

test_om = movmean(rad2deg(kep(:,5)), 500);
figure()
plot(test_om)
grid on
title('Argument of Periapsis Movmean Test');
xlabel('time [s]'); ylabel('\omega [°]');

test_theta = movmean(rad2deg(kep(:,6)), 500);
figure()
plot(test_theta)
grid on
title('True Anomaly Movmean Test');
xlabel('time [s]'); ylabel('\theta [°]');


%% All plots together

figure()
plot(a_vect)
grid on
hold on 
plot(kep(:,1))
plot(test_a)
hold off
title('Semi major axis');
legend('Cartesian','Gauss','Filtered');
xlabel('time [s]'); ylabel('a [km]');

figure()
plot(e_vect)
grid on
hold on 
plot(kep(:,2))
plot(test_e)
hold off
title('Eccentricity');
legend('Cartesian','Gauss','Filtered');
xlabel('time [s]'); ylabel('e [-]');

figure()
plot(rad2deg(i_vect))
grid on
hold on 
plot(rad2deg(kep(:,3)))
plot(test_i)
hold off
title('Inclination');
legend('Cartesian','Gauss','Filtered');
xlabel('time [s]'); ylabel('i [°]');

figure()
plot(rad2deg(Om_vect))
grid on
hold on
plot(rad2deg(kep(:,4)))
plot(test_Om)
plot(rad2deg(Om_sec))
hold off
title('RAAN');
legend('Cartesian','Gauss','Filtered','Secular');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(om_vect))
grid on
hold on 
plot(rad2deg(kep(:,5)))
plot(test_om)
plot(rad2deg(om_sec))
hold off
title('Argument of Periapsis');
legend('Cartesian','Gauss','Filtered','Secular');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(theta_vect))
grid on
hold on 
plot(rad2deg(kep(:,6)))
plot(test_theta)
plot(rad2deg(theta_sec))
hold off
title('True Anomaly');
legend('Cartesian','Gauss','Filtered','Secular');
xlabel('time [s]'); ylabel('\theta [°]');

%% ground track plot (perturbed)

om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)

orbit_number = orbit.ratio_k*10;
tspan_dim = 100000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart_perturbed(y0, tspan*orbit_number, earth.mu, theta_g, om_E, earth.r, earth.J2, spacecraft.AM, spacecraft.cD);
groundTrackPlot(lon, lat, "EarthTexture.jpg")

orbit.a_rep = aFinder(orbit.ratio_k, orbit.ratio_m, om_E, earth.mu);
T = 2*pi*sqrt( orbit.a_rep^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car(orbit.a_rep, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart_perturbed(y0, tspan*orbit_number, earth.mu, theta_g, om_E, earth.r, earth.J2, spacecraft.AM, spacecraft.cD);
groundTrackPlot(lon, lat, "EarthTexture.jpg")


%% other celestial body

addpath("Functions\")
addpath("Functions_custom\")

clear, clc
close all
A = importdata("test.csv");

% e_vect = A.data(:, 1);
% plot(movmean(e_vect, 50))
% grid on

% orbit data
orbit_new_object.a = A.data(1,10);
orbit_new_object.e = A.data(1,1);
orbit_new_object.i = deg2rad(A.data(1,3));
orbit_new_object.OM = deg2rad(A.data(1,4)); %Taken from similar object
orbit_new_object.om = deg2rad(A.data(1,5)); %Taken from similar object
orbit_new_object.theta=deg2rad(A.data(1,9)); %Taken from similar object
orbit_new_object.kep = [orbit_new_object.a orbit_new_object.e orbit_new_object.i orbit_new_object.OM orbit_new_object.om orbit_new_object.theta];
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

%% perturbations of new body


kep_body = [A.data(1,10), A.data(1,1), deg2rad(A.data(1,3)), deg2rad(A.data(1,4)), deg2rad(A.data(1,5)), deg2rad(A.data(1,9)), earth.mu];
[r0, v0] = kep2car(kep_body);

n_orbits = 30;
n_points = length(A.data);

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan= linspace( 0, T*n_orbits, n_points );

y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

%Plot orbit perturbed
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
plot(a_vect)
grid on
title('Semi major axis');
xlabel('time [s]'); ylabel('a [km]');

figure()
plot(e_vect)
grid on
title('Eccentricity');
xlabel('time [s]'); ylabel('e [-]');

figure()
plot(rad2deg(i_vect))
grid on
title('Inclination');
xlabel('time [s]'); ylabel('i [°]');

figure()
plot(rad2deg(Om_vect))
grid on
title('RAAN');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(om_vect))
grid on
title('Argument of periapsis');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(theta_vect))
grid on
title('True Anomaly');
xlabel('time [s]'); ylabel('\theta [°]');


% 
% figure
% plot(movmean(A.data(:,1), 50))
% hold on
% plot(movmean(e_vect, 65))

%% J2 Secular Effect

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );


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
    Om_sec = parameter*cos(i_vect(i)); 
    om_sec =parameter*(((5/2)*(sin(i_vect(i)))^2)-2); 
    theta_sec = parameter_2*(1-((2/3)*((sin(i_vect(i)))^2))); 

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


%% Perturbation over time

kep0_new_object = [orbit_new_object.a; orbit_new_object.e; orbit_new_object.i; orbit_new_object.OM; orbit_new_object.om; orbit_new_object.theta];
[t, kep_new_object] = ode113(@(t,kep_new_object) Pert_guass_eq_tnh_frame(t, kep_new_object, @(t,kep_new_object) acc_pert_fun(t, kep_new_object, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, kep0_new_object, options);

figure()
plot(kep_new_object(:,1))
grid on
title('Semi major axis');
xlabel('time [s]'); ylabel('a [km]');

figure()
plot(kep_new_object(:,2))
grid on
title('Eccentricity');
xlabel('time [s]'); ylabel('e [-]');

figure()
plot(rad2deg(kep_new_object(:,3)))
grid on
title('Inclination');
xlabel('time [s]'); ylabel('i [°]');

figure()
plot(rad2deg(kep_new_object(:,4)))
grid on
title('RAAN');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(kep_new_object(:,5)))
grid on
title('Argument of Periapsis');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(kep_new_object(:,6)))
grid on
title('True Anomaly');
xlabel('time [s]'); ylabel('\theta [°]');



%% test for movmean


test_a = movmean(kep_new_object(:,1), 500);
figure()
plot(test_a)
grid on
title('Semi major axis Movmean Test');

test_e= movmean(kep_new_object(:,2), 500);
figure()
plot(test_e)
grid on
title('Eccentricity Movmean Test');

test_i = movmean(rad2deg(kep_new_object(:,3)), 500);
figure()
plot(test_i)
grid on
title('Inclination Movmean Test');

test_Om = movmean(rad2deg(kep_new_object(:,4)), 500);
figure()
plot(test_Om)
grid on
title('RAAN Movmean Test');

test_om = movmean(rad2deg(kep_new_object(:,5)), 500);
figure()
plot(test_om)
grid on
title('Argument of Periapsis Movmean Test');

test_theta = movmean(rad2deg(kep_new_object(:,6)), 500);
figure()
plot(test_theta)
grid on
title('True Anomaly Movmean Test');



%% All plots together

figure()
plot(a_vect)
grid on
hold on 
plot(kep_new_object(:,1))
plot(test_a)
hold off
title('Semi major axis');
legend('Cartesian','Gauss','Filtered');
xlabel('time [s]'); ylabel('a [km]');

figure()
plot(e_vect)
grid on
hold on 
plot(kep_new_object(:,2))
plot(test_e)
hold off
title('Eccentricity');
legend('Cartesian','Gauss','Filtered');
xlabel('time [s]'); ylabel('e [-]');

figure()
plot(rad2deg(i_vect))
grid on
hold on 
plot(rad2deg(kep_new_object(:,3)))
plot(test_i)
hold off
title('Inclination');
legend('Cartesian','Gauss','Filtered');
xlabel('time [s]'); ylabel('i [°]');

figure()
plot(rad2deg(Om_vect))
grid on
hold on
plot(rad2deg(kep_new_object(:,4)))
plot(test_Om)
plot(rad2deg(Om_sec))
hold off
title('RAAN');
legend('Cartesian','Gauss','Filtered','Secular');
xlabel('time [s]'); ylabel('\Omega [°]');

figure()
plot(rad2deg(om_vect))
grid on
hold on 
plot(rad2deg(kep_new_object(:,5)))
plot(test_om)
plot(rad2deg(om_sec))
hold off
title('Argument of Periapsis');
legend('Cartesian','Gauss','Filtered','Secular');
xlabel('time [s]'); ylabel('\omega [°]');

figure()
plot(rad2deg(theta_vect))
grid on
hold on 
plot(rad2deg(kep_new_object(:,6)))
plot(test_theta)
plot(rad2deg(theta_sec))
hold off
title('True Anomaly');
legend('Cartesian','Gauss','Filtered','Secular');
xlabel('time [s]'); ylabel('\theta [°]');

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

