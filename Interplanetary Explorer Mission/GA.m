%% test

clear, clc
close all

addpath("Functions\")
addpath("Functions\time\")
addpath("Functions_custom\")
addpath("GA functions\")

time_instant = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.plnetId = 5;
arrival.bodyId = 79;

impact_parameter = 100;
plane_hyp = [0 0 1];


time_instant_mjd200 = date2mjd2000(time_instant);

[departure.kep, mu_sun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.plnetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);


% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% orbit points
[departure.r0, departure.v0] = kep2car([departure.kep, mu_sun]);
[flyby.r0, flyby.v0] = kep2car([flyby.kep, mu_sun]);
[arrival.r0, arrival.v0] = kep2car([arrival.kep, mu_sun]);

% hyp plane, test
hyp_plane = plane3points(departure.r0, flyby.r0);
flyby.r0_new = newpoint(flyby.r0, departure.r0, impact_parameter, hyp_plane);







%% newpoint check

v_impact = [1 8 0];
ref_point = [3 2 0];
impact_param = 2;
hyp_planeVect = [0 0 1];

[v_new] = newpoint(v_impact, ref_point, impact_param, hyp_planeVect);

plot3(v_impact(1), v_impact(2), v_impact(3), 'or')
hold on
plot3(ref_point(1), ref_point(2), ref_point(3), 'or')
plot3(v_new(1), v_new(2), v_new(3), 'xr')
axis equal, grid on

check = dot((v_new - ref_point), (v_impact - v_new));
impact_check = norm(v_impact - v_new);

