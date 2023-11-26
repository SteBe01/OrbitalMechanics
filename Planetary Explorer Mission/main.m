%% orbit data

addpath("Functions\")

clear, clc
close all

% orbit data - goup 2346
orbit.a = 0.8016 * 1e4;
orbit.e = 0.1678;
orbit.i = deg2rad(50.3442);
orbit.ratio_k = 12;
orbit.ratio_m = 1;

orbit.OM = deg2rad(0);
orbit.om = deg2rad(0);

% perturbation: J2 and Drag (cD = 2.1, A/M = 0.0171 m^2/kg)


%% unperturbed 2bp orbit

mu = astroConstants(13);

T = 2*pi*sqrt( orbit.a^3/mu );
tspan= linspace( 0, T, 100 );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, 0, mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp(t,y,mu), tspan, y0, options );

figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

title('Two-body problem orbit');
axis equal, grid on,hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or' )


%% ground track plot (unperturbed)

R_E = astroConstants(23);
om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)
mu = astroConstants(13);
orbit_number = orbit.ratio_k*10;
tspan_dim = 100000;

T = 2*pi*sqrt( orbit.a^3/mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, 0, mu);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, mu, theta_g, om_E);
groundTrackPlot(lon, lat, "EarthTexture.jpg")

orbit.a_rep = aFinder(orbit.ratio_k, orbit.ratio_m, om_E, mu);
T = 2*pi*sqrt( orbit.a_rep^3/mu );
tspan= linspace( 0, T, tspan_dim );
[r0, v0] = kep2car(orbit.a_rep, orbit.e, orbit.i, orbit.OM, orbit.om, 0, mu);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, mu, theta_g, om_E);
groundTrackPlot(lon, lat, "EarthTexture.jpg")


%% perturbations

mu = astroConstants(13);
Re = astroConstants(23);
om_E = deg2rad(15.04) / 3600;
A_M = 0.0171;
J2 = astroConstants(9);
cD = 2.1;

T = 2*pi*sqrt( orbit.a^3/mu );
tspan= linspace( 0, T*23, 10000 );

[r0, v0] = kep2car(orbit.a, orbit.e, orbit.i, orbit.OM, orbit.om, 0, mu);
y0 = [ r0'; v0' ];

options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
[ T, Y ] = ode113( @(t,y) ode_2bp_perturbed( t, y, mu, Re, J2, om_E, A_M, cD), tspan, y0, options );

figure()
plot3( Y(:,1), Y(:,2), Y(:,3), '-' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

title('Two-body problem orbit');
axis equal, grid on, hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or' )
plot3( Y(end,1), Y(end,2), Y(end,3), 'or' )

a_vect = [];
e_vect = [];
for i = 1:length(Y)
    [a, e, ~, ~, ~, ~] = car2kep(Y(i,1:3), Y(i,4:6), mu);
    a_vect = [a_vect a];
    e_vect = [e_vect e];
end

figure
plot(a_vect)
grid on
figure
plot(e_vect)
grid on

