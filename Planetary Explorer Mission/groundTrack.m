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
plot3( Y(:,1), Y(:,2), Y(:,3),'r')
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');

title('Two-body problem orbit');
axis equal, grid on,hold on
earthPlot;
plot3( Y(1,1), Y(1,2), Y(1,3), 'or','MarkerEdgeColor','blue' )


%% ground track plot (unperturbed)

om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)

orbit_number = orbit.ratio_k;
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
[r0, v0] = kep2car(orbit.a_rep, orbit.e, orbit.i, orbit.OM, orbit.om, orbit.theta, earth.mu);
y0 = [ r0'; v0' ];
[~, ~, lon, lat] = groundTrack_cart(y0, tspan*orbit_number, earth.mu, theta_g, om_E);
groundTrackPlot(lon, lat, "EarthTexture.jpg")


%% ground track plot (perturbed)

om_E = 15.04;                   % deg/h
theta_g = 0;                    % theta Greenwich (t0)

orbit_number = orbit.ratio_k;
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

