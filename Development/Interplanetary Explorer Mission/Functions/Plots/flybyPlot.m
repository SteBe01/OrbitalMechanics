function [data, image, lgd] = flybyPlot(departure_time, flyby_time, arrival_time, departure_Id, flyby_Id, arrival_Id, flybyTime)

% Creates the plot for the flyby
%
% Usage
% [data, image, lgd] = flybyPlot(departure_time, flyby_time, arrival_time, departure_Id, flyby_Id, arrival_Id, flybyTime)
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
% data                  [-]     all data available      [-]
% image                 [1x1]   handler for the image   [-]
% lgd                   [1x1]   handler for the legend  [-]
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

orbitType = 0;

ToF_dep_flyby = (flyby_time - departure_time) * 24 * 60 * 60;
ToF_flyby_arr = (arrival_time-flyby_time) * 24 * 60 * 60;

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
mu_planet = astroConstants(flyby_Id + 10);

if flyby_Id <= 9
    m_planet = astroConstants(flyby_Id + 10) / astroConstants(1);
    m_sun = astroConstants(4) / astroConstants(1);
else
    error("r_orb_planet not defined!")
end

[departure.r0, ~] = kep2car([departure.kep, ksun]);
[flyby.r0, flyby.v0] = kep2car([flyby.kep, ksun]);
[arrival.r0, ~] = kep2car([arrival.kep, ksun]);

r_soi = norm(flyby.r0) * (m_planet/m_sun)^(2/5);

% first transfer ARC
[A_1, ~, ~, ERROR1, ~, VF_1, ~, ~] = lambertMR(departure.r0, flyby.r0, ToF_dep_flyby, ksun, orbitType, 0);
% second transfer ARC
[A_2, ~, ~, ERROR2, VI_2, ~, ~, ~] = lambertMR(flyby.r0, arrival.r0, ToF_flyby_arr, ksun, orbitType, 0);

if A_1 < 0 || A_2 < 0
    error("A_1 < 0 || A_2 < 0")
end
if ERROR1 || ERROR2
    error("Error inside lambert")
end

% Heliocentric velocities
v_inf_minus = VF_1;
v_inf_plus = VI_2;

% planetocentric velocities
V_inf_minus = v_inf_minus - flyby.v0;
V_inf_plus = v_inf_plus - flyby.v0;

[rp, flag] = rpsolver(V_inf_minus, V_inf_plus, flyby_Id);
if flag
    error("rpsolver returned an error")
end

a_minus = -mu_planet/norm(V_inf_minus)^2;
a_plus = -mu_planet/norm(V_inf_plus)^2;

e_minus = 1 - rp/a_minus;
e_plus = 1 - rp/a_plus;
delta_minus = 2*asin(1/e_minus);
delta_plus = 2*asin(1/e_plus);

delta = atan2((norm(cross(V_inf_minus, V_inf_plus)))/(norm(V_inf_minus)*norm(V_inf_plus)), (dot(V_inf_minus, V_inf_plus))/(norm(V_inf_minus)*norm(V_inf_plus)));
if abs((delta_plus + delta_minus)/2 - delta) > 0.0001
    error("(delta_plus + delta_minus)/2 ~= delta")
end

impact_param_minus = -a_minus/tan(delta/2);
impact_param_plus = -a_plus/tan(delta/2);

hyp_plane = cross(V_inf_minus, V_inf_plus);
hyp_plane = hyp_plane ./ norm(hyp_plane);

rp_direction = rodrigues_rotation(V_inf_minus./norm(V_inf_minus), delta_minus/2 - pi/2, hyp_plane);
rp_vect = rp .* rp_direction;

Vp_hyp_minus = sqrt(norm(V_inf_minus)^2+2*mu_planet/rp);
Vp_hyp_plus = sqrt(norm(V_inf_plus)^2+2*mu_planet/rp);
deltaVp = Vp_hyp_plus - Vp_hyp_minus;

deltaVtot = norm(v_inf_plus - v_inf_minus);

directionPeri = rodrigues_rotation(rp_direction, pi/2, hyp_plane);

% options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

y0 = [rp_vect directionPeri.*Vp_hyp_minus];
tspan = linspace( 0, -flybyTime/2, 1000 );
[ ~, Y_minus ] = ode113( @(t,y) ode_2bp(t,y,mu_planet), tspan, y0, options );
Y_len = length(Y_minus);

[Y_minus, maxIndex] = vectorCut(Y_minus, r_soi);
if maxIndex < Y_len
    tspan = linspace( 0, tspan(maxIndex), 1000 );
    [ ~, Y_minus ] = ode113( @(t,y) ode_2bp(t,y,mu_planet), tspan, y0, options );
end

image = figure();
plot3( Y_minus(:,1), Y_minus(:,2), Y_minus(:,3), LineWidth=2)
hold on, grid on, axis equal

y0 = [rp_vect directionPeri.*Vp_hyp_plus];
tspan = linspace( 0, flybyTime/2, 1000 );
[ ~, Y_plus ] = ode113( @(t,y) ode_2bp(t,y,mu_planet), tspan, y0, options );

[Y_plus, maxIndex] = vectorCut(Y_plus, r_soi);
if maxIndex < Y_len
    tspan = linspace( 0, tspan(maxIndex), 1000 );
    [ ~, Y_plus ] = ode113( @(t,y) ode_2bp(t,y,mu_planet), tspan, y0, options );
end

plot3( Y_plus(:,1), Y_plus(:,2), Y_plus(:,3), LineWidth=2)

% asymptotes
V_inf_minus_direction = V_inf_minus./norm(V_inf_minus);
V_inf_plus_direction = V_inf_plus./norm(V_inf_plus);

asyStart1 = (rp-a_minus) .* rp_direction;
mult1 = vecnorm(Y_minus(end, 1:3) - asyStart1);
plot3(asyStart1(1), asyStart1(2), asyStart1(3), 'xr', HandleVisibility='off')
asy_minus = [-V_inf_minus_direction*mult1 + asyStart1; asyStart1];
plot3(asy_minus(:,1), asy_minus(:,2), asy_minus(:,3), '--b', LineWidth=1)

asyStart2 = (rp-a_plus) .* rp_direction;
mult2 = vecnorm(Y_plus(end, 1:3) - asyStart2);
plot3(asyStart2(1), asyStart2(2), asyStart2(3), 'xr', HandleVisibility='off')
asy_plus = [V_inf_plus_direction*mult2 + asyStart2; asyStart2];
plot3(asy_plus(:,1), asy_plus(:,2), asy_plus(:,3), '--r', LineWidth=1)

if norm(asy_minus(2,:)) < norm(asy_plus(2,:))
    point = asy_plus(2,:);
else
    point = asy_minus(2,:);
end

point = [point; -point];
plot3(point(:,1), point(:,2), point(:,3), '--m', LineWidth=1)

% planet velocity
% planetV_plot = [flyby.v0; 0 0 0].*1000000;
% plot3(planetV_plot(:,1), planetV_plot(:,2), planetV_plot(:,3), '--m', LineWidth=1)

if hyp_plane(3) > 0
    view(hyp_plane)
else
    view(-hyp_plane)
end

[x1, y1, z1] = sphere;
mult = astroConstants(flyby_Id + 20);
surface(x1*mult, y1*mult, z1*mult, EdgeColor="none", FaceColor=[0 0 0]);

PlanetNames = ["Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto", "Sun"];

% SoI plot
if vecnorm(Y_plus(end, 1:3)) > r_soi || vecnorm(Y_minus(end, 1:3)) > r_soi
    [x1, y1, z1] = sphere;
    surf(x1 * r_soi, y1 * r_soi, z1 * r_soi, 'FaceAlpha', .2, 'EdgeColor', 'none', FaceColor=[0 1 0]);

    lgd = legend("Hyperbola -Inf", "Hyperbola +Inf", "Asyptote -Inf", "Asyptote +Inf", "Hyperbola axis", PlanetNames(flyby_Id), "Sphere of Influence");
else
    lgd = legend("Hyperbola -Inf", "Hyperbola +Inf", "Asyptote -Inf", "Asyptote +Inf", "Hyperbola axis", PlanetNames(flyby_Id));
end
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title("Flyby hyperbola")

% time of flight sphere of influence
deltaT1 = tof_soi(a_plus, e_plus, r_soi, mu_planet);
deltaT2 = tof_soi(a_minus, e_minus, r_soi, mu_planet);

data.rp = rp;
data.r_planet = mult;
data.h = rp - mult;
data.deltaTot = (delta_minus + delta_plus)/2;
data.rsoi = r_soi;
data.deltaVp = deltaVp;
data.deltaVtot = deltaVtot;
data.ratio = deltaVp/deltaVtot;
data.ToF = deltaT1 + deltaT2;

data.a.plus = a_plus;
data.a.minus = a_minus;
data.e.plus = e_plus;
data.e.minus = e_minus;
data.delta.plus = delta_plus;
data.delta.minus = delta_minus;
data.impact_param.plus = impact_param_plus;
data.impact_param.minus = impact_param_minus;
data.tof.plus = deltaT1;
data.tof.minus = deltaT2;

end


%% functions

function [vect, i] = vectorCut(vect, maxVal)
    for i = 1:length(vect)
        if vecnorm(vect(i, 1:3)) > maxVal
            break
        end
    end
    if i <= 5
        error("reduce flybyTime")
    end
    vect = vect(1:i, :);
end

