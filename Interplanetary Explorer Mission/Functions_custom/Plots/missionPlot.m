function [image, lgd] = missionPlot(departure_time, flyby_time, arrival_time, departure_Id, flyby_Id, arrival_Id)

% Creates the plot for the whole mission
%
% Usage
% [image, lgd] = missionPlot(departure_time, flyby_time, arrival_time, departure_Id, flyby_Id, arrival_Id)
%
% Input arguments:
% ----------------------------------------------------------------
% departure_time        [1x1]   departure time          [date]
% flyby_time            [1x1]   flyby time              [date]
% arrival_time          [1x1]   arrival time            [date]
% departure_Id          [1x1]   Id of dep. body         [-]
% flyby_Id              [1x1]   Id of flyby body        [-]
% arrival_Id            [1x1]   Id of arr. body         [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% image                 [1x1]   handler for the image   [-]
% lgd                   [1x1]   handler for the legend  [-]

orbitType = 0;

ToF_dep_flyby = (flyby_time - departure_time) * 24 * 60 * 60;
ToF_flyby_arr = (arrival_time-flyby_time) * 24 * 60 * 60;

if ToF_dep_flyby <= 0 || ToF_flyby_arr <= 0
    error("ToF is not valid!")
end

if departure_Id <= 10
    [departure.kep, ~] = uplanet(departure_time, departure_Id);
else
    [departure.kep, ~] = ephNEO(departure_time, departure_Id);
end
if flyby_Id <= 10
    [flyby.kep, ~] = uplanet(flyby_time, flyby_Id);
else
    [flyby.kep, ~] = ephNEO(flyby_time, flyby_Id);
end
if arrival_Id <= 10
    [arrival.kep, ~] = uplanet(arrival_time, arrival_Id);
else
    [arrival.kep, ~] = ephNEO(arrival_time, arrival_Id);
end

ksun = astroConstants(4);
PlanetNames = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Sun"];

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% orbit propagation - departure
[departure.r0, departure.v0] = kep2car([departure.kep, ksun]);
departure.y0 = [departure.r0 departure.v0];

departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [1/s]
departure.tspan= linspace( 0, departure.T_orb, 200 );
[ ~, departure.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), departure.tspan, departure.y0, options );

% orbit propagation - flyby
[flyby.r0, flyby.v0] = kep2car([flyby.kep, ksun]);
flyby.y0 = [flyby.r0 flyby.v0];

flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [1/s]
flyby.tspan= linspace( 0, flyby.T_orb, 200 );
[ ~, flyby.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), flyby.tspan, flyby.y0, options );

% orbit propagation - arrival
[arrival.r0, arrival.v0] = kep2car([arrival.kep, ksun]);
arrival.y0 = [arrival.r0 arrival.v0];

arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [1/s]
arrival.tspan= linspace( 0, arrival.T_orb, 200 );
[ ~, arrival.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), arrival.tspan, arrival.y0, options );

% first transfer ARC
[A_1, ~, ~, ERROR_1, VI_1, ~, ~, ~] = lambertMR(departure.r0, flyby.r0, ToF_dep_flyby, ksun, orbitType, 0);
y0_1 = [ departure.r0 VI_1 ];
% set time span
tspan_1 = linspace( 0,ToF_dep_flyby, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% perform the integration
[ ~, Y_1 ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_1, y0_1, options );

% second transfer ARC
[A_2, ~, ~, ERROR_2, VI_2, ~, ~, ~] = lambertMR(flyby.r0, arrival.r0, ToF_flyby_arr, ksun, orbitType, 0);
y0_2 = [ flyby.r0 VI_2 ];
% set time span
tspan_2 = linspace( 0,ToF_flyby_arr, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% perform the integration
[ ~, Y_2 ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan_2, y0_2, options );

if ERROR_1 || ERROR_2
    error("Lambert returned an error")
end
if A_1 < 0 || A_2 < 0
    error("A_1 < 0 || A_2 < 0")
end

% legend settings
if departure_Id <= 10
    departure_name = PlanetNames(departure_Id);
else
    departure_name = strcat("Asteroid N.", num2str(departure_Id));
end
if flyby_Id <= 10
    flyby_name = PlanetNames(flyby_Id);
else
    flyby_name = strcat("Asteroid N.", num2str(flyby_Id));
end
if arrival_Id <= 10
    arrival_name = PlanetNames(arrival_Id);
else
    arrival_name = strcat("Asteroid N.", num2str(arrival_Id));
end

% plot the results
image = figure();
plot3(departure.Y(:,1), departure.Y(:,2), departure.Y(:,3), '-', 'color', 'b', LineWidth=1)
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Complete mission');
axis equal, grid on, hold on
plot3(flyby.Y(:,1),flyby.Y(:,2), flyby.Y(:,3), '-', 'color', 'r', LineWidth=1)
plot3(arrival.Y(:,1), arrival.Y(:,2), arrival.Y(:,3), '-', 'color', 'g', LineWidth=1)
plot3(Y_1(:,1), Y_1(:,2), Y_1(:,3), '--', 'color', 'm', 'HandleVisibility', 'off', LineWidth=1)
plot3(Y_2(:,1), Y_2(:,2), Y_2(:,3), '--', 'color', 'm', 'HandleVisibility', 'off', LineWidth=1)
plot3(departure.r0(1), departure.r0(2), departure.r0(3), 'o', 'Color', 'b', 'MarkerFaceColor', 'b')
plot3(flyby.r0(1), flyby.r0(2), flyby.r0(3), 'o', 'Color', 'r', 'MarkerFaceColor', 'r')
plot3(arrival.r0(1), arrival.r0(2), arrival.r0(3), 'o', 'Color', 'g', 'MarkerFaceColor', 'g')
plot3(0, 0, 0, 'o', 'Color', 'y', 'MarkerFaceColor', 'y')
lgd = legend(strcat(departure_name, ' orbit'), strcat(flyby_name, ' orbit'), strcat(arrival_name, ' orbit'), departure_name, flyby_name, arrival_name, 'Sun');
hold off

end

