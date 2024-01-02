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

n_orbits = orbit.ratio_k*8;
n_points = 10000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan = linspace( 0, T*n_orbits, n_points );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

figure()
plot3( Y(:,1), Y(:,2), Y(:,3),'r')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag');
axis equal, grid on, hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or',MarkerEdgeColor='b' )
plot3( Y(end,1), Y(end,2), Y(end,3), 'or',MarkerEdgeColor='b' )


a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
i_vect = zeros(length(Y), 1);
Om_vect = zeros(length(Y), 1);
om_vect = zeros(length(Y), 1);
theta_vect = zeros(length(Y), 1);

for i = 1:length(Y)
    [a_vect(i), e_vect(i), i_vect(i), Om_vect(i), om_vect(i), theta_vect(i)] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
end

Period=2*pi*sqrt( orbit.a^3/earth.mu );

% figure()
% plot(tspan./(Period), a_vect,'b')
% grid on
% title('Semi major axis Evolution');
% xlabel('time [T]'); ylabel('a [km]');
% 
% figure()
% plot(tspan./(Period), e_vect,'b')
% grid on
% title('Eccentricity Evolution');
% xlabel('time [T]'); ylabel('e [-]');
% 
% figure()
% plot(tspan./(Period), rad2deg(i_vect),'b')
% grid on
% title('Inclination Evolution');
% xlabel('time [T]'); ylabel('i [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(Om_vect),'b')
% grid on
% title('RAAN');
% xlabel('time [T]'); ylabel('\Omega [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(om_vect),'b')
% grid on
% title('Argument of Periapsis Evolution');
% xlabel('time [T]'); ylabel('\omega [°]');

% figure()
% plot(tspan./(Period), rad2deg(unwrap(theta_vect)),'b')
% grid on
% title('True Anomaly Evolution');
% xlabel('time [T]'); ylabel('\theta [°]');


%% perturbations - Gauss's planetary equations

s0 = [orbit.a; orbit.e; orbit.i; orbit.OM; orbit.om; orbit.theta];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, kep] = ode113(@(t,s) eq_motion(t, s, @(t,s) acc_pert_fun(t, s, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, s0, options);

Period=2*pi*sqrt( orbit.a^3/earth.mu );


figure()
plot(tspan./(Period), kep(:,1))
grid on
title('Semi major axis Evolution');
xlabel('time [T]'); ylabel('a [km]');

figure()
plot(tspan./(Period), kep(:,2))
grid on
title('Eccentricity Evolution');
xlabel('time [T]'); ylabel('e [-]');

figure()
plot(tspan./(Period), rad2deg(kep(:,3)))
grid on
title('Inclination Evolution');
xlabel('time [T]'); ylabel('i [°]');

figure()
plot(tspan./(Period), rad2deg(kep(:,4)))
grid on
title('RAAN Evolution');
xlabel('time [T]'); ylabel('\Omega [°]');

figure()
plot(tspan./(Period), rad2deg(kep(:,5)))
grid on
title('Argument of Periapsis Evolution');
xlabel('time [T]'); ylabel('\omega [°]');

figure()
plot(tspan./(Period), kep(:,6))
grid on
title('True Anomaly Evolution');
xlabel('time [T]'); ylabel('\theta [°]');


%% filter for movmean

n_orbits = orbit.ratio_k*8;

num_elements=length(kep(:,1));
num_elements_per_period=num_elements/n_orbits;

filter_a = movmean(kep(:,1), num_elements_per_period);
filter_e= movmean(kep(:,2), num_elements_per_period);
filter_i = movmean(rad2deg(kep(:,3)), num_elements_per_period);
filter_Om = movmean(rad2deg(kep(:,4)), num_elements_per_period);
filter_om = movmean(rad2deg(kep(:,5)), num_elements_per_period);
filter_theta = movmean(rad2deg(kep(:,6)), num_elements_per_period);


% figure()
% plot(filter_a,'b')
% grid on
% title('Semi major axis Movmean filter');
% xlabel('time [s]'); ylabel('a [km]');
% 
% figure()
% plot(filter_e)
% grid on
% title('Eccentricity Movmean filter');
% xlabel('time [s]'); ylabel('e [-]');
% 
% figure()
% plot(filter_i)
% grid on
% title('Inclination Movmean filter');
% xlabel('time [s]'); ylabel('i [°]');
% 
% figure()
% plot(filter_Om)
% grid on
% title('RAAN Movmean filter');
% xlabel('time [s]'); ylabel('\Omega [°]')
% 
% figure()
% plot(filter_om)
% grid on
% title('Argument of Periapsis Movmean filter');
% xlabel('time [s]'); ylabel('\omega [°]');

% figure()
% plot(filter_theta)
% grid on
% title('True Anomaly Movmean filter');
% xlabel('time [s]'); ylabel('\theta [°]');


%% All plots together

%% Semi major axis

Period=2*pi*sqrt( orbit.a^3/earth.mu );

% Only filtered
figure()
plot(tspan./(Period),movmean(a_vect,num_elements_per_period),'m')
grid on
hold on 
plot(tspan./(Period),filter_a,'r')
grid on
hold off
title('Semi major axis Evolution Filtered');
legend('Cartesian','Gauss','Location', 'Best');
xlabel('time [T]'); ylabel('a [km]');

% Filtered vs unfiltered
figure()
plot(tspan./(Period),a_vect,'m')
grid on
hold on 
plot(tspan./(Period),kep(:,1),'b')
plot(tspan./(Period),movmean(a_vect,num_elements_per_period),'y')
plot(tspan./(Period),filter_a,'r')
grid on
hold off
title('Semi major axis Evolution Filtered vs Unfltered');
legend('Cartesian','Gauss','Cartesian Filtered','Gauss Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('a [km]');


%% Eccentricity

% Only filtered
figure()
plot(tspan./(Period),movmean(e_vect,num_elements_per_period),'m')
grid on
hold on 
plot(tspan./(Period),filter_e,'r')
grid on
hold off
title('Eccentricity Evolution Filtered');
legend('Cartesian','Gauss','Location', 'Best');
xlabel('time [T]'); ylabel('e [-]');

%  Filtered vs unfiltered
figure()
plot(tspan./(Period),e_vect,'m')
grid on
hold on 
plot(tspan./(Period),kep(:,2),'b')
plot(tspan./(Period),movmean(e_vect,num_elements_per_period),'y')
plot(tspan./(Period),filter_e,'r')
grid on
hold off
title('Eccentricity Evolution');
legend('Cartesian','Gauss','Cartesian Filtered','Gauss Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('e [-]');

%% Inclination 

% Only filtered
figure()
plot(tspan./(Period),rad2deg(movmean(i_vect,num_elements_per_period)),'m')
grid on
hold on 
plot(tspan./(Period),filter_i,'r')
grid on
hold off
title('Inclination Evolution');
legend('Cartesian','Gauss','Location', 'Best');
xlabel('time [T]'); ylabel('i [°]');


% Filtered vs unfiltered 
figure()
plot(tspan./(Period),rad2deg(wrapTo2Pi(i_vect)),'m')
grid on
hold on 
plot(tspan./(Period),rad2deg(wrapTo2Pi(kep(:,3))),'b')
plot(tspan./(Period),rad2deg(wrapTo2Pi(movmean(i_vect,num_elements_per_period))),'y')
plot(tspan./(Period),filter_i,'r')
grid on
hold off
title('Inclination Evolution');
legend('Cartesian','Gauss','Cartesian Filtered','Gauss Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('i [°]');

%% RAAN

% Only filtered
figure()
plot(tspan./(Period),rad2deg(movmean(Om_vect,num_elements_per_period)),'m')
grid on
hold on
plot(tspan./(Period),filter_Om,'r')
grid on
hold off
title('RAAN Evolution');
legend('Cartesian','Gauss','Location', 'Best');
xlabel('time [T]'); ylabel('\Omega [°]');

% Filtered vs Unfiltered
figure()
plot(tspan./(Period),rad2deg(Om_vect),'m')
grid on
hold on
plot(tspan./(Period),rad2deg(kep(:,4)),'b')
plot(tspan./(Period),rad2deg(movmean(Om_vect,num_elements_per_period)),'y')
plot(tspan./(Period),filter_Om,'r')
grid on
hold off
title('RAAN Evolution');
legend('Cartesian','Gauss','Cartesian Filtered','Gauss Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('\Omega [°]');

%% Argument of Periapsis 

% Only filtered
figure()
plot(tspan./(Period),rad2deg(movmean(om_vect,num_elements_per_period)),'m')
grid on
hold on 
plot(tspan./(Period),filter_om,'r')
grid on
hold off
title('Argument of Periapsis Evolution');
legend('Cartesian','Gauss','Location', 'Best');
xlabel('time [T]'); ylabel('\omega [°]');

% Filtered vs unfiltered
figure()
plot(tspan./(Period),rad2deg(om_vect),'m')
grid on
hold on 
plot(tspan./(Period),rad2deg(kep(:,5)),'b')
plot(tspan./(Period),rad2deg(movmean(om_vect,num_elements_per_period)),'y')
plot(tspan./(Period),filter_om,'r')
grid on
hold off
title('Argument of Periapsis Evolution');
legend('Cartesian','Gauss','Cartesian Filtered','Gauss Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('\omega [°]');

%% True Anomaly

% Only filtered
figure()
plot(tspan./(Period),rad2deg(movmean(unwrap(theta_vect),num_elements_per_period)),'m')
grid on
hold on 
plot(tspan./(Period),rad2deg(filter_theta),'r')
grid on
hold off
title('True Anomaly Evolution');
legend('Cartesian','Gauss','Location', 'Best');
xlabel('time [T]'); ylabel('\theta [°]');


figure()
plot(tspan./(Period),rad2deg(unwrap(theta_vect)),'b')
grid on
hold on 
plot(tspan./(Period),rad2deg(kep(:,6)),'b')
plot(tspan./(Period),rad2deg(movmean(unwrap(theta_vect),num_elements_per_period)),'m')
plot(tspan./(Period),rad2deg(filter_theta),'r')
grid on
hold off
title('True Anomaly Evolution');
legend('Cartesian','Gauss','Cartesian Filtered','Gauss Filtered','Location', 'Best');
xlabel('time [T]'); ylabel('\theta [°]');



%% Difference Certesian and Gauss

figure()
a_diff=abs(a_vect-kep(:,1))/orbit.a;
plot(tspan./(Period),a_diff,'b')
grid on
title('Semi major axis Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|a_c_a_r_t - a_g_a_u_s_s|/ a_0 [km]');

figure()
e_diff=abs(e_vect-kep(:,2));
plot(tspan./(Period),e_diff,'b')
grid on
title('Eccentricity Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|e_c_a_r_t - e_g_a_u_s_s| [-]');

figure()
i_diff=abs(i_vect-kep(:,3))/(2*pi());
plot(tspan./(Period),rad2deg(wrapTo2Pi(i_diff)),'b')
grid on
title('Inclination Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|i_c_a_r_t - i_g_a_u_s_s|/2 \pi [°]');

figure()
Om_diff=abs(Om_vect-kep(:,4))/(2*pi());
plot(tspan./(Period),rad2deg(wrapTo2Pi(Om_diff)),'b')
grid on
title('RAAN Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|\Omega_c_a_r_t - \Omega_g_a_u_s_s|/2 \pi [°]');

figure()
om_diff=abs(om_vect-kep(:,5))/(2*pi());
plot(tspan./(Period),rad2deg(wrapTo2Pi(om_diff)),'b')
grid on
title('Argument of Periapsis Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|\omega_c_a_r_t - \omega_g_a_u_s_s|/2 \pi [°]');


figure()
theta_diff=abs(unwrap(theta_vect)-kep(:,6))./abs(kep(:,6));
plot(tspan./(Period),rad2deg(wrapTo2Pi(theta_diff)),'b')
grid on
title('True Anomaly Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|\theta_c_a_r_t - \theta_g_a_u_s_s|/ \theta_g_a_u_s_s[°]');


%% Simulation 

n_orbits = orbit.ratio_k*8;
n_points = 100;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan = linspace( 0, T*n_orbits, n_points );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), tspan, y0, options );

figure()
plot3( Y(:,1), Y(:,2), Y(:,3),'r')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag');
axis equal, grid on, hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or',MarkerEdgeColor='b' )
plot3( Y(end,1), Y(end,2), Y(end,3), 'or',MarkerEdgeColor='b' )

figure()
plot3( Y(1,1), Y(1,2), Y(1,3), 'or',MarkerEdgeColor='b' )
hold on
plot3( Y(end,1), Y(end,2), Y(end,3), 'or',MarkerEdgeColor='b' )
ani1=animatedline('Color','r');
for i=1:length(tspan)
    addpoints(ani1,Y(i,1), Y(i,2), Y(i,3));
    earthPlot;
    drawnow
    pause(0.01);
end

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


%% without propagation

orbit_new_object.a_no_prop = A.data(:,10);
orbit_new_object.e_no_prop = A.data(:,1);
orbit_new_object.i_no_prop = A.data(:,3);
orbit_new_object.OM_no_prop = A.data(:,4);
orbit_new_object.om_no_prop = A.data(:,5);
orbit_new_object.theta_no_prop = A.data(:,9);
 
n_orbits = orbit.ratio_k*8;
n_points = length(A.data);

T = 2*pi*sqrt( orbit_new_object.a^3/earth.mu );
tspan= linspace( 0, T*n_orbits, n_points );
Period=2*pi*sqrt( orbit_new_object.a^3/earth.mu );

% figure()
% plot(tspan./(Period), orbit_new_object.a_no_prop)
% grid on
% title('Semi major axis');
% xlabel('time [T]'); ylabel('a [km]');
% 
% figure()
% plot(tspan./(Period), orbit_new_object.e_no_prop)
% grid on
% title('Eccentricity');
% xlabel('time [T]'); ylabel('e [-]');
% 
% figure()
% plot(tspan./(Period), orbit_new_object.i_no_prop)
% grid on
% title('Inclination');
% xlabel('time [T]'); ylabel('i [°]');
% 
% figure()
% plot(tspan./(Period), orbit_new_object.OM_no_prop)
% grid on
% title('RAAN');
% xlabel('time [T]'); ylabel('\Omega [°]');
% 
% figure()
% plot(tspan./(Period), orbit_new_object.om_no_prop)
% grid on
% title('Argument of periapsis');
% xlabel('time [T]'); ylabel('\omega [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(unwrap(orbit_new_object.theta_no_prop)))
% grid on
% title('True Anomaly');
% xlabel('time [T]'); ylabel('\theta [°]');


%% perturbations - cartesian coordinates

kep_body = [A.data(1,10), A.data(1,1), deg2rad(A.data(1,3)), deg2rad(A.data(1,4)), deg2rad(A.data(1,5)), deg2rad(A.data(1,9)), earth.mu];
[r0, v0] = kep2car(kep_body);

n_orbits = orbit.ratio_k*8;
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

% figure()
% plot(tspan./(Period), a_vect)
% grid on
% title('Semi major axis');
% xlabel('time [T]'); ylabel('a [km]');
% 
% figure()
% plot(tspan./(Period), e_vect)
% grid on
% title('Eccentricity');
% xlabel('time [T]'); ylabel('e [-]');
% 
% figure()
% plot(tspan./(Period), rad2deg(i_vect))
% grid on
% title('Inclination');
% xlabel('time [T]'); ylabel('i [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(Om_vect))
% grid on
% title('RAAN');
% xlabel('time [T]'); ylabel('\Omega [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(om_vect))
% grid on
% title('Argument of periapsis');
% xlabel('time [T]'); ylabel('\omega [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(unwrap(theta_vect)))
% grid on
% title('True Anomaly');
% xlabel('time [T]'); ylabel('\theta [°]');


%% perturbations - Gauss's planetary equations

s0 = [orbit_new_object.a; orbit_new_object.e; orbit_new_object.i; orbit_new_object.OM; orbit_new_object.om; orbit_new_object.theta];
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[t, kep_new_object] = ode113(@(t,s) eq_motion(t, s, @(t,s) acc_pert_fun(t, s, earth.mu, earth.r, earth.J2, earth.om, spacecraft.AM, spacecraft.cD), earth.mu), tspan, s0, options);

% figure()
% plot(tspan./(Period), kep_new_object(:,1),'b')
% grid on
% title('Semi major axis');
% xlabel('time [T]'); ylabel('a [km]');
% 
% figure()
% plot(tspan./(Period), kep_new_object(:,2),'b')
% grid on
% title('Eccentricity');
% xlabel('time [T]'); ylabel('e [-]');
% 
% figure()
% plot(tspan./(Period), rad2deg(kep_new_object(:,3)),'b')
% grid on
% title('Inclination');
% xlabel('time [T]'); ylabel('i [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(kep_new_object(:,4)),'b')
% grid on
% title('RAAN');
% xlabel('time [T]'); ylabel('\Omega [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(kep_new_object(:,5)),'b')
% grid on
% title('Argument of Periapsis');
% xlabel('time [T]'); ylabel('\omega [°]');
% 
% figure()
% plot(tspan./(Period), rad2deg(unwrap(kep_new_object(:,6))),'b')
% grid on
% title('True Anomaly');
% xlabel('time [T]'); ylabel('\theta [°]');


%% filter for movmean

n_orbits = orbit.ratio_k*8;

num_elements=length(kep_new_object(:,1));
num_elements_per_period=num_elements/n_orbits;

filter_a = movmean(kep_new_object(:,1), num_elements_per_period);
filter_e= movmean(kep_new_object(:,2), num_elements_per_period);
filter_i = movmean(rad2deg(kep_new_object(:,3)), num_elements_per_period);
filter_Om = movmean(rad2deg(kep_new_object(:,4)), num_elements_per_period);
filter_om = movmean(rad2deg(kep_new_object(:,5)), num_elements_per_period);
filter_theta =  movmean(rad2deg(unwrap(kep_new_object(:,6))), num_elements_per_period);


% figure()
% plot(filter_a,'b')
% grid on
% title('Semi major axis Movmean filter');
% xlabel('time [s]'); ylabel('a [km]');
% 
% figure()
% plot(filter_e)
% grid on
% title('Eccentricity Movmean filter');
% xlabel('time [s]'); ylabel('e [-]');
% 
% figure()
% plot(filter_i)
% grid on
% title('Inclination Movmean filter');
% xlabel('time [s]'); ylabel('i [°]');
% 
% figure()
% plot(filter_Om)
% grid on
% title('RAAN Movmean filter');
% xlabel('time [s]'); ylabel('\Omega [°]')
% 
% figure()
% plot(filter_om)
% grid on
% title('Argument of Periapsis Movmean filter');
% xlabel('time [s]'); ylabel('\omega [°]');
% 
% figure()
% plot(filter_theta)
% grid on
% title('True Anomaly Movmean filter');
% xlabel('time [s]'); ylabel('\theta [°]');

%% All plots together

%% Semi major axis

Period=2*pi*sqrt( orbit_new_object.a^3/earth.mu );

figure()
plot(tspan./(Period),orbit_new_object.a_no_prop,'g')
grid on
hold on
plot(tspan./(Period),a_vect,'c')
plot(tspan./(Period),kep_new_object(:,1),'b')
plot(tspan./(Period), movmean(orbit_new_object.a_no_prop,num_elements_per_period),'m')
plot(tspan./(Period),movmean(a_vect,num_elements_per_period),'y')
plot(tspan./(Period),filter_a,'r')
hold off
title('Semi major axis Evolution Real Body');
legend('Real Data','Cartesian','Gauss', ...
    'Real Data Filtered','Cartesian Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('a [km]');



%% Eccentricity

figure()
plot(tspan./(Period),orbit_new_object.e_no_prop,'g')
grid on
hold on
plot(tspan./(Period),e_vect,'c')
plot(tspan./(Period),kep_new_object(:,2),'b')
plot(tspan./(Period), movmean(orbit_new_object.e_no_prop,num_elements_per_period),'m')
plot(tspan./(Period),movmean(e_vect,num_elements_per_period),'y')
plot(tspan./(Period),filter_e,'r')
hold off
title('Eccentricity Evolution Real Body');
legend('Real Data','Cartesian','Gauss', ...
    'Real Data Filtered','Cartesian Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('e [-]');


%% Inclination 

figure()
plot(tspan./(Period),orbit_new_object.i_no_prop,'g')
grid on
hold on
plot(tspan./(Period),rad2deg(i_vect),'c')
plot(tspan./(Period),rad2deg(kep_new_object(:,3)),'b')
plot(tspan./(Period), movmean(orbit_new_object.i_no_prop,num_elements_per_period),'m')
plot(tspan./(Period),rad2deg(movmean(i_vect,num_elements_per_period)),'y')
plot(tspan./(Period),filter_i,'r')
hold off
title('Inclination Evolution Real Body');
legend('Real Data','Cartesian','Gauss', ...
    'Real Data Filtered','Cartesian Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('i [°]');

%% RAAN

figure()
plot(tspan./(Period),unwrap(orbit_new_object.OM_no_prop),'g')
grid on
hold on
plot(tspan./(Period),rad2deg(Om_vect),'c')
plot(tspan./(Period),rad2deg(kep_new_object(:,4)),'b')
plot(tspan./(Period), movmean(unwrap(orbit_new_object.OM_no_prop),num_elements_per_period),'m')
plot(tspan./(Period),rad2deg(movmean(Om_vect,num_elements_per_period)),'y')
plot(tspan./(Period),filter_Om,'r')
hold off
title('RAAN Evolution Real Body');
legend('Real Data','Cartesian','Gauss', ...
    'Real Data Filtered','Cartesian Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('\Omega [°]');

%% Argument of Periapsis 

figure()
plot(tspan./(Period),unwrap(orbit_new_object.om_no_prop),'g')
grid on
hold on
plot(tspan./(Period),rad2deg(om_vect),'c')
plot(tspan./(Period),rad2deg(kep_new_object(:,5)),'b')
plot(tspan./(Period), movmean(unwrap(orbit_new_object.om_no_prop),num_elements_per_period),'m')
plot(tspan./(Period),rad2deg(movmean(om_vect,num_elements_per_period)),'y')
plot(tspan./(Period),filter_om,'r')
hold off
title('Argument of periapsis Evolution Real Body');
legend('Real Data','Cartesian','Gauss', ...
    'Real Data Filtered','Cartesian Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('\omega [°]');

%% True Anomaly

figure()
plot(tspan./(Period),rad2deg(unwrap(orbit_new_object.theta_no_prop)),'g')
grid on
hold on
plot(tspan./(Period),rad2deg(unwrap(theta_vect)),'c')
plot(tspan./(Period),rad2deg(kep_new_object(:,6)),'b')
plot(tspan./(Period), rad2deg(movmean(unwrap(orbit_new_object.theta_no_prop),num_elements_per_period)),'m')
plot(tspan./(Period),rad2deg(movmean(unwrap(theta_vect),num_elements_per_period)),'y')
plot(tspan./(Period),filter_theta,'r')
hold off
title('True Anomaly Evolution Real Body');
legend('Real Data','Cartesian','Gauss', ...
    'Real Data Filtered','Cartesian Filtered','Gauss Filtered', ...
    'Location','Best');
xlabel('time [T]'); ylabel('\theta [°]');


%% Difference Real Data and Gauss

figure()
a_diff=abs(orbit_new_object.a_no_prop-kep_new_object(:,1))/orbit_new_object.a;
plot(tspan./(Period),a_diff,'b')
grid on
title('Semi major axis Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|a_r_e_a_l - a_g_a_u_s_s|/ a_0 [km]');

figure()
e_diff=abs(orbit_new_object.e_no_prop-kep_new_object(:,2));
plot(tspan./(Period),e_diff,'b')
grid on
title('Eccentricity Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|e_r_e_a_l - e_g_a_u_s_s| [-]');

figure()
i_diff=abs(deg2rad(orbit_new_object.i_no_prop)-kep_new_object(:,3))/(2*pi());
plot(tspan./(Period),rad2deg(i_diff),'b')
grid on
title('Inclination Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|i_r_e_a_l - i_g_a_u_s_s|/2 \pi [°]');

figure()
Om_diff=abs(deg2rad(unwrap(orbit_new_object.OM_no_prop))-unwrap(kep_new_object(:,4)))/(2*pi());
plot(tspan./(Period),rad2deg(Om_diff),'b')
grid on
title('RAAN Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|\Omega_r_e_a_l - \Omega_g_a_u_s_s|/2 \pi [°]');

figure()
om_diff=abs(deg2rad(unwrap(orbit_new_object.om_no_prop))-unwrap(kep_new_object(:,5)))/(2*pi());
plot(tspan./(Period),rad2deg(om_diff),'b')
grid on
title('Argument of Periapsis Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|\omega_r_e_a_l - \omega_g_a_u_s_s|/2 \pi [°]');


figure()
theta_diff=abs(unwrap(orbit_new_object.theta_no_prop)-unwrap(kep_new_object(:,6)))./abs(unwrap(orbit_new_object.theta_no_prop));
plot(tspan./(Period),rad2deg(wrapTo2Pi(theta_diff)),'b')
grid on
title('True Anomaly Evolution Propagation Method Difference');
xlabel('time [T]'); ylabel('|\theta_r_e_a_l - \theta_g_a_u_s_s|/ \theta_r_e_a_l[°]');


