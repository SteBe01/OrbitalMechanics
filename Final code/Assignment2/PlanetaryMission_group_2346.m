%% data - group 2346

% This script was used to generate all the plots for the Planetary Mission
%
% SECTIONS:
% line 54  - Unperturbed 2bp orbit
% line 79  - Ground track unperturbed (36 orbits) + unperturbed (1 orbit)
% line 106 - Ground track perturbed (36 orbits) + unperturbed (36 orbits)
% line 126 - Ground track repetition (new a, 36 orbits), with perturbed reprtition (new a, 36 orbits)
% line 152 - Perturbations - Cartesian coordinates
% line 184 - Perturbations - Cartesian's planetary equations
% line 199 - Perturbations - Gauss's planetary equations
% line 206 - Filter for movmean
% line 219 - Evolution of Keplerian Elements (Gaussian method) and filtered
% line 294 - Difference Certesian and Gauss
% line 370 - Real celestial body
% line 414 - Plot Perturbation New Body
% line 445 - Perturbations - Cartesian coordinates
% line 461 - Perturbations - Cartesian's planetary equations
% line 475 - Perturbations - Gauss's planetary equations
% line 482 - Filter for movmean
% line 495 - All plots together
% line 590 - Differences Real Data and Gauss
% line 666 - Test animation 

restoredefaultpath
addpath(genpath("."))

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


%% Unperturbed 2bp orbit

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, 100 );

[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,earth.mu), tspan, y0, options );

% Orbit plot
figure()
plot3( Y(1,1), Y(1,2), Y(1,3), 'or','MarkerEdgeColor','blue' ,LineWidth=2, MarkerSize=7)
hold on
plot3( Y(:,1), Y(:,2), Y(:,3),'r',LineWidth=2)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit');
axis equal, grid on,hold on
s = earthPlot;
legend("Start", Location="best")
rotate(s, [0 0 1],-70,[0 0 0]);
view(-40,20)


%% Ground track unperturbed (36 orbits) + unperturbed (1 orbit)

om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)

% 36 Orbits
orbit_number = orbit.ratio_k*3;
tspan_dim = 100000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot(lon, lat, "EarthTexture.jpg")
title('Ground Track - 36 Orbits and 1 Orbit Without Perturbations');

% 1 Orbit
orbit_number = 1;
T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.kep, earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot2(lon, lat, "red", 3)


%% Ground track perturbed (36 orbits) + unperturbed (36 orbits)

% 36 Orbits
orbit_number = orbit.ratio_k*3;
tspan_dim = 100000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.kep earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart_perturbed(y0, tspan*orbit_number, earth.mu, theta_g, om_E, earth.r, earth.J2, spacecraft.AM, spacecraft.cD);
groundTrackPlot(lon, lat, "EarthTexture.jpg")
title('Ground Track - 36 Orbits With and Without Perturbations');

[r0, v0] = kep2car([orbit.kep earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot2(lon, lat, "red", 1.5)


%% Ground track repetition (new a, 36 orbits), with perturbed reprtition (new a, 36 orbits)

% 36 Orbits
orbit_number = orbit.ratio_k*3;
tspan_dim = 100000;

orbit.a_rep = aFinder(orbit.ratio_k, orbit.ratio_m, om_E, earth.mu);

T = 2*pi*sqrt( orbit.a_rep^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.a_rep orbit.kep(2:end) earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart_perturbed(y0, tspan*orbit_number, earth.mu, theta_g, om_E, earth.r, earth.J2, spacecraft.AM, spacecraft.cD);
groundTrackPlot(lon, lat, "EarthTexture.jpg")
title('Repeating Ground Track - 36 Orbits With and Without Perturbations');

[r0, v0] = kep2car([orbit.a_rep orbit.kep(2:end) earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot2(lon, lat, "red", 1.5)

disp("The semi-major axis given was 8016 km")
fprintf("The semi-major axis needed for a repeating ground track with a ratio of 12:1 is: " + orbit.a_rep + " km\n")


%% Perturbations - Cartesian coordinates

n_orbits = orbit.ratio_k*5;
n_points = 100000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan = linspace( 0, T*n_orbits, n_points );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

% Orbit Propagation
figure()
plot3( Y(1,1), Y(1,2), Y(1,3), 'or','MarkerEdgeColor','blue' ,LineWidth=2, MarkerSize=7)
hold on
plot3( Y(end,1), Y(end,2), Y(end,3), 'or','MarkerEdgeColor','red' ,LineWidth=2, MarkerSize=7)
scatter3( Y(:,1), Y(:,2), Y(:,3),1,1:length(tspan))
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag');
axis equal, grid on, hold on
s = earthPlot;
hold off

rotate(s, [0 0 1],-70,[0 0 0]);
view(-40,20)

legend("Start", "End", Location="best")


%% Perturbations - Cartesian's planetary equations

Period = 2*pi*sqrt(orbit.a^3/earth.mu);
a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
i_vect = zeros(length(Y), 1);
Om_vect = zeros(length(Y), 1);
om_vect = zeros(length(Y), 1);
theta_vect = zeros(length(Y), 1);

for i = 1:length(Y)
    [a_vect(i), e_vect(i), i_vect(i), Om_vect(i), om_vect(i), theta_vect(i)] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
end


%% Perturbations - Gauss's planetary equations

s0 = [orbit.a; orbit.e; orbit.i; orbit.OM; orbit.om; orbit.theta];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, kep] = ode113(@(t,s) eq_motion(t, s, @(t,s) acc_pert_fun(t, s, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, s0, options);


%% Filter for movmean

num_elements=length(kep(:,1));
num_elements_per_period=num_elements/n_orbits;

filter_a = movmean(kep(:,1), num_elements_per_period,'Endpoints','fill');
filter_e= movmean(kep(:,2), num_elements_per_period ,'Endpoints','fill');
filter_i = movmean(rad2deg(kep(:,3)), num_elements_per_period,'Endpoints','fill');
filter_Om = movmean(rad2deg(kep(:,4)), num_elements_per_period,'Endpoints','fill');
filter_om = movmean(rad2deg(kep(:,5)), num_elements_per_period,'Endpoints','fill');
filter_theta = movmean(rad2deg(unwrap(kep(:,6))), num_elements_per_period,'Endpoints','fill');


%% Evolution of Keplerian Elements (Gaussian method) and filtered

% Semi major axis
figure()
grid on
hold on 
plot(tspan./(Period),kep(:,1),'b')
plot(tspan./(Period),filter_a,'r',LineWidth=2)
grid on
hold off
title('Semi major axis Evolution Filtered vs Unfltered');
legend('Gaussian','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('a [km]');

% Eccentricity
figure()
grid on
hold on 
plot(tspan./(Period),kep(:,2),'b')
plot(tspan./(Period),filter_e,'r',LineWidth=2)
grid on
hold off
title('Eccentricity Evolution');
legend('Gaussian','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('e [-]');

% Inclination 
figure()
grid on
hold on 
plot(tspan./(Period),rad2deg(wrapTo2Pi(kep(:,3))),'b')
plot(tspan./(Period),filter_i,'r',LineWidth=2)
grid on
hold off
title('Inclination Evolution');
legend('Gaussian','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('i [°]');

% RAAN
figure()
grid on
hold on
plot(tspan./(Period),rad2deg(kep(:,4)),'b',LineWidth=1.5)
plot(tspan./(Period),filter_Om,'r',LineWidth=1)
grid on
hold off
title('RAAN Evolution');
legend('Gaussian','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('\Omega [°]');

% Argument of Periapsis 
figure()
grid on
hold on 
plot(tspan./(Period),rad2deg(kep(:,5)),'b',LineWidth=1.5)
plot(tspan./(Period),filter_om,'r',LineWidth=2)
grid on
hold off
title('Argument of Periapsis Evolution');
legend('Gaussian','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('\omega [°]');

% True Anomaly
figure()
grid on
hold on 
plot(tspan./(Period),rad2deg(kep(:,6)),'b',LineWidth=1.5)
plot(tspan./(Period),rad2deg(movmean(unwrap(theta_vect),num_elements_per_period)),'r',LineWidth=1)
grid on
hold off
title('True Anomaly Evolution');
legend('Gaussian','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('\theta [°]');


%% Difference Certesian and Gauss

% Semi-major axis
figure()
a_diff=abs(a_vect-kep(:,1))/orbit.a;
plot(tspan./(Period),a_diff,'b')
hold on
plot(tspan./(Period),movmean(a_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Semi major axis Evolution Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|a_c_a_r_t - a_g_a_u_s_s|/ a_0 [-]');

% Eccentricity
figure()
e_diff=abs(e_vect-kep(:,2));
plot(tspan./(Period),e_diff,'b')
hold on
plot(tspan./(Period),movmean(e_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Eccentricity Evolution Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|e_c_a_r_t - e_g_a_u_s_s| [-]');


% Inclination
figure()
i_diff=abs(i_vect-kep(:,3))/(2*pi());
plot(tspan./(Period),i_diff,'b')
hold on
plot(tspan./(Period),movmean(i_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Inclination Evolution Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|i_c_a_r_t - i_g_a_u_s_s|/2 \pi [-]');

% RAAN
figure()
Om_diff=abs(Om_vect-kep(:,4))/(2*pi());
plot(tspan./(Period),Om_diff,'b')
hold on
plot(tspan./(Period),movmean(Om_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('RAAN Evolution Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|\Omega_c_a_r_t - \Omega_g_a_u_s_s|/2 \pi [-]');

% Argument of periapsis
figure()
om_diff=abs(om_vect-kep(:,5))/(2*pi());
plot(tspan./(Period),om_diff,'b')
hold on
plot(tspan./(Period),movmean(om_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Argument of Periapsis Evolution Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|\omega_c_a_r_t - \omega_g_a_u_s_s|/2 \pi [-]');

% True Anomaly
figure()
theta_diff=abs(unwrap(theta_vect)-kep(:,6))./abs(kep(:,6));
plot(tspan./(Period),theta_diff,'b')
hold on
plot(tspan./(Period),movmean(theta_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('True Anomaly Evolution Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|\theta_c_a_r_t - \theta_g_a_u_s_s|/ \theta_g_a_u_s_s[-]');


%% Real celestial body

clear

A = importdata("EXPRESS-MD2.csv");

% orbit data
orbit_new_object.a = A.data(1,10);
orbit_new_object.e = A.data(1,1);
orbit_new_object.i = deg2rad(A.data(1,3));
orbit_new_object.OM = deg2rad(A.data(1,4)); 
orbit_new_object.om = deg2rad(A.data(1,5)); 
orbit_new_object.theta = deg2rad(A.data(1,9)); 
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

% without propagation
n_orbits = orbit.ratio_k*5;

orbit_new_object.a_no_prop = A.data(:,10);
orbit_new_object.e_no_prop = A.data(:,1);
orbit_new_object.i_no_prop = A.data(:,3);
orbit_new_object.OM_no_prop = A.data(:,4);
orbit_new_object.om_no_prop = A.data(:,5);
orbit_new_object.theta_no_prop = A.data(:,9);

n_points = length(A.data);

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan_nb= linspace( 0, T*n_orbits, n_points );
Period=2*pi*sqrt( orbit_new_object.a^3/earth.mu );


%% Plot Perturbation New Body

kep_body = [A.data(1,10), A.data(1,1), deg2rad(A.data(1,3)), deg2rad(A.data(1,4)), deg2rad(A.data(1,5)), deg2rad(A.data(1,9)), earth.mu];
[r0, v0] = kep2car(kep_body);

n_points = 100000;

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan_nb= linspace( 0, T*n_orbits, n_points );

y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan_nb, y0, options );

% Orbit Propagation
figure()
plot3( Y(1,1), Y(1,2), Y(1,3), 'or','MarkerEdgeColor','blue' ,LineWidth=2, MarkerSize=7)
hold on
plot3( Y(end,1), Y(end,2), Y(end,3), 'or','MarkerEdgeColor','red' ,LineWidth=2, MarkerSize=7)
scatter3( Y(:,1), Y(:,2), Y(:,3),1,1:length(tspan_nb))
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag New Body');
axis equal, grid on, hold on
s = earthPlot;
hold off
rotate(s, [0 0 1],-70,[0 0 0]);
view(-40,20)
legend("Start", "End", Location="best")


%% Perturbations - Cartesian coordinates

kep_body = [A.data(1,10), A.data(1,1), deg2rad(A.data(1,3)), deg2rad(A.data(1,4)), deg2rad(A.data(1,5)), deg2rad(A.data(1,9)), earth.mu];
[r0, v0] = kep2car(kep_body);

n_points = length(A.data);

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan_nb= linspace( 0, T*n_orbits, n_points );

y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan_nb, y0, options );


%% Perturbations - Cartesian's planetary equations

a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
i_vect = zeros(length(Y), 1);
Om_vect = zeros(length(Y), 1);
om_vect = zeros(length(Y), 1);
theta_vect = zeros(length(Y), 1);

for i = 1:length(Y)
    [a_vect(i), e_vect(i), i_vect(i), Om_vect(i), om_vect(i), theta_vect(i)] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
end


%% Perturbations - Gauss's planetary equations

s0 = [orbit_new_object.a; orbit_new_object.e; orbit_new_object.i; orbit_new_object.OM; orbit_new_object.om; orbit_new_object.theta];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, kep_new_object] = ode113(@(t,s) eq_motion(t, s, @(t,s) acc_pert_fun(t, s, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan_nb, s0, options);


%% filter for movmean

num_elements=length(kep_new_object(:,1));
num_elements_per_period=num_elements/n_orbits;

filter_a = movmean(kep_new_object(:,1), num_elements_per_period,'Endpoints','fill');
filter_e= movmean(kep_new_object(:,2), num_elements_per_period,'Endpoints','fill');
filter_i = movmean(rad2deg(kep_new_object(:,3)), num_elements_per_period,'Endpoints','fill');
filter_Om = movmean(rad2deg(kep_new_object(:,4)), num_elements_per_period,'Endpoints','fill');
filter_om = movmean(rad2deg(kep_new_object(:,5)), num_elements_per_period,'Endpoints','fill');
filter_theta =  movmean(rad2deg(unwrap(kep_new_object(:,6))), num_elements_per_period,'Endpoints','fill');


%% All plots together

% Semi major axis
Period=2*pi*sqrt( orbit_new_object.a^3/earth.mu );

figure()
plot(tspan_nb./(Period),orbit_new_object.a_no_prop,'m')
grid on
hold on
plot(tspan_nb./(Period),kep_new_object(:,1),'b')
plot(tspan_nb./(Period), movmean(orbit_new_object.a_no_prop,num_elements_per_period,'Endpoints','fill'),'g',LineWidth=1.5)
plot(tspan_nb./(Period),filter_a,'r',LineWidth=1.5)
hold off
title('Semi major axis Evolution Real Body');
legend('Real Data','Gauss', ...
    'Real Data Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('a [km]');

% Eccentricity
figure()
plot(tspan_nb./(Period),orbit_new_object.e_no_prop,'m')
grid on
hold on
plot(tspan_nb./(Period),kep_new_object(:,2),'b')
plot(tspan_nb./(Period), movmean(orbit_new_object.e_no_prop,num_elements_per_period,'Endpoints','fill'),'g',LineWidth=1.5)
plot(tspan_nb./(Period),filter_e,'r',LineWidth=1.5)
hold off
title('Eccentricity Evolution Real Body');
legend('Real Data','Gauss', ...
    'Real Data Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('e [-]');

% Inclination
figure()
plot(tspan_nb./(Period),orbit_new_object.i_no_prop,'m')
grid on
hold on
plot(tspan_nb./(Period),rad2deg(kep_new_object(:,3)),'b')
plot(tspan_nb./(Period), movmean(orbit_new_object.i_no_prop,num_elements_per_period,'Endpoints','fill'),'g',LineWidth=1.5)
plot(tspan_nb./(Period),filter_i,'r',LineWidth=1.5)
hold off
title('Inclination Evolution Real Body');
legend('Real Data','Gauss', ...
    'Real Data Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('i [°]');

% RAAN
figure()
plot(tspan_nb./(Period),unwrap(orbit_new_object.OM_no_prop),'m',LineWidth=1.2)
grid on
hold on
plot(tspan_nb./(Period),rad2deg(kep_new_object(:,4)),'b',LineWidth=1.2)
plot(tspan_nb./(Period), movmean(unwrap(orbit_new_object.OM_no_prop),num_elements_per_period,'Endpoints','fill'),'g')
plot(tspan_nb./(Period),filter_Om,'r')
hold off
title('RAAN Evolution Real Body');
legend('Real Data','Gauss', ...
    'Real Data Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('\Omega [°]');

% Argument of Periapsis 
figure()
plot(tspan_nb./(Period),unwrap(orbit_new_object.om_no_prop),'m',LineWidth=1.2)
grid on
hold on
plot(tspan_nb./(Period),rad2deg(kep_new_object(:,5)),'b',LineWidth=1.2)
plot(tspan_nb./(Period), movmean(unwrap(orbit_new_object.om_no_prop),num_elements_per_period,'Endpoints','fill'),'g')
plot(tspan_nb./(Period),filter_om,'r')
hold off
title('Argument of periapsis Evolution Real Body');
legend('Real Data','Gauss', ...
    'Real Data Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('\omega [°]');

% True Anomaly
figure()
plot(tspan_nb./(Period),rad2deg(unwrap(orbit_new_object.theta_no_prop)),'m')
grid on
hold on
plot(tspan_nb./(Period),rad2deg(kep_new_object(:,6)),'b',LineWidth=1.2)
plot(tspan_nb./(Period), rad2deg(movmean(unwrap(orbit_new_object.theta_no_prop),num_elements_per_period,'Endpoints','fill')),'g')
plot(tspan_nb./(Period),filter_theta,'r')
hold off
title('True Anomaly Evolution Real Body');
legend('Real Data','Gauss', ...
    'Real Data Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('\theta [°]');


%% Differences Real Data and Gauss

% Semi-major axis
figure()
a_diff=abs(orbit_new_object.a_no_prop-kep_new_object(:,1))/orbit_new_object.a;
plot(tspan_nb./(Period),a_diff,'b')
hold on
plot(tspan_nb./(Period), movmean(a_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Semi major axis Evolution Real Data vs  Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|a_r_e_a_l - a_g_a_u_s_s|/ a_0 [-]');

% Eccentricity
figure()
e_diff=abs(orbit_new_object.e_no_prop-kep_new_object(:,2));
plot(tspan_nb./(Period),e_diff,'b')
hold on
plot(tspan_nb./(Period), movmean(e_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Eccentricity Evolution Real Data vs Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|e_r_e_a_l - e_g_a_u_s_s| [-]');


% Inclination
figure()
i_diff=abs(deg2rad(orbit_new_object.i_no_prop)-kep_new_object(:,3))/(2*pi());
plot(tspan_nb./(Period),i_diff,'b')
hold on
plot(tspan_nb./(Period), movmean(i_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('Inclination Evolution Real Data vs Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|i_r_e_a_l - i_g_a_u_s_s|/2 \pi [-]');

% RAAN
figure()
Om_diff=abs(deg2rad(unwrap(orbit_new_object.OM_no_prop))-unwrap(kep_new_object(:,4)))/(2*pi());
plot(tspan_nb./(Period),Om_diff,'b',LineWidth=1.2)
hold on
plot(tspan_nb./(Period), movmean(Om_diff,num_elements_per_period,'Endpoints','fill'),'r')
hold off
grid on
title('RAAN Evolution Evolution Real Data vs Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|\Omega_r_e_a_l - \Omega_g_a_u_s_s|/2 \pi [-]');

% Argument of periapsis
figure()
om_diff=abs(deg2rad(unwrap(orbit_new_object.om_no_prop))-unwrap(kep_new_object(:,5)))/(2*pi());
plot(tspan_nb./(Period),om_diff,'b',LineWidth=1)
hold on
plot(tspan_nb./(Period), movmean(om_diff,num_elements_per_period,'Endpoints','fill'),'r')
hold off
grid on
title('Argument of Periapsis Evolution Real Data vs Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|\omega_r_e_a_l - \omega_g_a_u_s_s|/2 \pi [-]');

% True anomaly
figure()
theta_diff=abs(unwrap(orbit_new_object.theta_no_prop)-unwrap(kep_new_object(:,6)))./abs(unwrap(orbit_new_object.theta_no_prop));
plot(tspan_nb./(Period),theta_diff,'b')
hold on
plot(tspan_nb./(Period), movmean(theta_diff,num_elements_per_period,'Endpoints','fill'),'r',LineWidth=1.5)
hold off
grid on
title('True Anomaly Evolution Real Data vs Propagation Method Difference');
legend('Difference','Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('|\theta_r_e_a_l - \theta_g_a_u_s_s|/ \theta_r_e_a_l[-]');


%% Test animation 

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

n_orbits = orbit.ratio_k*5;
n_points = 10000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan = linspace( 0, T*n_orbits, n_points );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

%Orbit Plot
figure()
plot3( Y(1,1), Y(1,2), Y(1,3), 'or',MarkerEdgeColor='b' )
hold on
plot3( Y(end,1), Y(end,2), Y(end,3), 'or',MarkerEdgeColor='b' )
ani1=animatedline('Color','r');
earthPlot;
title('Animation')
for i=1:length(tspan)
    addpoints(ani1,Y(i,1), Y(i,2), Y(i,3));
    drawnow
end
