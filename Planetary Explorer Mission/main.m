%% data - group 2346

addpath("Functions\")

clear, clc
close all

% orbit data
orbit.a = 0.8016 * 1e4;
orbit.e = 0.1678;
orbit.i = deg2rad(50.3442);
orbit.OM = deg2rad(0);
orbit.om = deg2rad(0);
orbit.kep = [orbit.a orbit.e orbit.i orbit.OM orbit.om];
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

[r0, v0] = kep2car([orbit.kep, 0, earth.mu]);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ ~, Y ] = ode113( @(t,y) ode_2bp(t,y,earth.mu), tspan, y0, options );

figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
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
[r0, v0] = kep2car([orbit.kep, 0, earth.mu]);
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
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Two-body problem orbit, with J2 and air drag');
axis equal, grid on, hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or' )
plot3( Y(end,1), Y(end,2), Y(end,3), 'or' )

a_vect = zeros(length(Y), 1);
e_vect = zeros(length(Y), 1);
for i = 1:length(Y)
    [a_vect(i), e_vect(i), ~, ~, ~, ~] = car2kep(Y(i,1:3), Y(i,4:6), earth.mu);
end

figure
plot(a_vect)
grid on
title('a');
xlabel('time [s]'); ylabel('a [km]');

figure
plot(e_vect)
grid on
title('e');
xlabel('time [s]'); ylabel('e [-]');


%% test for movmean

test = movmean(a_vect, 500);

figure
plot(test)
grid on


%% ground track plot (perturbed)

om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)

orbit_number = orbit.ratio_k*10;
tspan_dim = 100000;

T = 2*pi*sqrt( orbit.a^3/earth.mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car([orbit.kep, 0, earth.mu]);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart_perturbed(y0, tspan*orbit_number, earth.mu, theta_g, om_E, earth.r, earth.J2, spacecraft.AM, spacecraft.cD);
groundTrackPlot(lon, lat, "EarthTexture.jpg")
