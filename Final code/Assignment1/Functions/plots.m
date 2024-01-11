%% plots

% This script contains code used to create plots.
% Run the "data" section (this one) first and then run each section
% independently.

clc, clear
close all

restoredefaultpath
addpath(genpath("."))

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058


%% orbits plot

clc
close all

time_instant = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(time_instant);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.planetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% orbit propagation - departure
[departure.r0, departure.v0] = kep2car([departure.kep, ksun]);
departure.y0 = [departure.r0 departure.v0];

departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [1/s]
tspan= linspace( 0, departure.T_orb, 200 );
[ ~, departure.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, departure.y0, options );

% orbit propagation - flyby
[flyby.r0, flyby.v0] = kep2car([flyby.kep, ksun]);
flyby.y0 = [flyby.r0 flyby.v0];

flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [1/s]
tspan= linspace( 0, flyby.T_orb, 200 );
[ ~, flyby.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, flyby.y0, options );

% orbit propagation - arrival
[arrival.r0, arrival.v0] = kep2car([arrival.kep, ksun]);
arrival.y0 = [arrival.r0 arrival.v0];

arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [1/s]
tspan= linspace( 0, arrival.T_orb, 200 );
[ ~, arrival.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, arrival.y0, options );

% plots
figure
axis equal, grid on, hold on
plot3(departure.Y(:, 1), departure.Y(:, 2), departure.Y(:, 3), Color="blu", HandleVisibility="off")
plot3(flyby.Y(:, 1), flyby.Y(:, 2), flyby.Y(:, 3), Color="green", HandleVisibility="off")
plot3(arrival.Y(:, 1), arrival.Y(:, 2), arrival.Y(:, 3), Color="red", HandleVisibility="off")

scatter3(departure.r0(1), departure.r0(2), departure.r0(3), "blu", "filled")
scatter3(flyby.r0(1), flyby.r0(2), flyby.r0(3), "green", "filled")
scatter3(arrival.r0(1), arrival.r0(2), arrival.r0(3), "red", "filled")
scatter3(0, 0, 0, "yellow", "filled")

legend("Departure (Saturn)", "Fly-By (Jupiter)", "Arrival (Asteroid N.79)", "Sun", Location="best")

xlim([min(min(min(departure.Y(:,1)), min(flyby.Y(:,1))), min(arrival.Y(:,1))) max(max(max(departure.Y(:,1)), max(flyby.Y(:,1))), max(arrival.Y(:,1)))])
ylim([min(min(min(departure.Y(:,2)), min(flyby.Y(:,2))), min(arrival.Y(:,2))) max(max(max(departure.Y(:,2)), max(flyby.Y(:,2))), max(arrival.Y(:,2)))])
zlim([min(min(min(departure.Y(:,3)), min(flyby.Y(:,3))), min(arrival.Y(:,3))) max(max(max(departure.Y(:,3)), max(flyby.Y(:,3))), max(arrival.Y(:,3)))])

xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")

view(30,30)


%% mission plot

clc, clear
close all

load("dvmin.mat")

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

[image, lgd] = missionPlot(xcust(1), xcust(2), xcust(3), departure.planetId, flyby.planetId, arrival.bodyId);
fontsize(image, 15, "points");
fontsize(lgd, 10, "points");


%% flyby plot

clc, clear
close all

load("dvmin.mat")
flybyTime = 1e7;

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

[data, image, lgd]= flybyPlot(xcust(1), xcust(2), xcust(3), departure.planetId, flyby.planetId, arrival.bodyId, flybyTime);
fontsize(image, 15, "points");
fontsize(lgd, 12, "points");

