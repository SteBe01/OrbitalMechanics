%% data

clc, clear
close all

addpath("Functions\")
addpath("Functions\time\")
addpath("Functions_custom\")

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058


%% orbit plot

clc
close all


time_instant = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.plnetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(time_instant);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.plnetId);
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


legend("Departure (Saturn)", "Fly-By (Jupiter)", "Arrival (Asteroid N.79)", Location="best")

xlim([min(min(min(departure.Y(:,1)), min(flyby.Y(:,1))), min(arrival.Y(:,1))) max(max(max(departure.Y(:,1)), max(flyby.Y(:,1))), max(arrival.Y(:,1)))])
ylim([min(min(min(departure.Y(:,2)), min(flyby.Y(:,2))), min(arrival.Y(:,2))) max(max(max(departure.Y(:,2)), max(flyby.Y(:,2))), max(arrival.Y(:,2)))])
zlim([min(min(min(departure.Y(:,3)), min(flyby.Y(:,3))), min(arrival.Y(:,3))) max(max(max(departure.Y(:,3)), max(flyby.Y(:,3))), max(arrival.Y(:,3)))])

xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")

view(30,30)


%% Find Departure and Arrival Dates Ranges

% ---Synodic Periods-----------

date_departure = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(date_departure);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.planetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);


departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [1/s]

flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [1/s]

arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [1/s]

%Synodic Period between departure planet and flyby
T_syn_dep_flyby=(departure.T_orb*flyby.T_orb)/(abs(departure.T_orb-flyby.T_orb)); %[s]
T_syn_dep_flyby_years=T_syn_dep_flyby/(60*60*24*astroConstants(32)); %[years]

%Synodic Period between flyby planet and arrival
T_syn_flyby_arr=(arrival.T_orb*flyby.T_orb)/(abs(arrival.T_orb-flyby.T_orb)); %[s]
T_syn_flyby_arr_years=T_syn_flyby_arr/(60*60*24*astroConstants(32)); %[years]

%Synodic Period between departure planet and arrival
T_syn_dep_arr=(arrival.T_orb*departure.T_orb)/(abs(arrival.T_orb-departure.T_orb)); %[s]
T_syn_dep_arr_years=T_syn_dep_arr/(60*60*24*astroConstants(32)); %[years]


%Set possible dates for departure and arrival
jd_arrival_flyby_planet=mjd20002jd(time_instant_mjd200)+T_syn_dep_flyby/(60*60*24);
date_flyby = jd2date(jd_arrival_flyby_planet);

jd_arrival_arrival_planet=jd_arrival_flyby_planet+T_syn_flyby_arr/(60*60*24);
date_arrival= jd2date(jd_arrival_arrival_planet);


%%
% ---------GA on Jupiter

%Heliocentric Velocities
V_minus=VF_1_actual+flyby.v0_actual; %[km/s]
V_plus=VI_2_actual+flyby.v0_actual; %[km/s]

%Planet Velocity
V_Planet=flyby.v0_actual;

%Planetocentric velocites
v_minus_inf=VF_1_actual;
v_plus_inf=VI_2_actual;


%Turning Angle
delta = acos(dot(v_plus_inf,v_minus_inf)/(norm(v_plus_inf)*norm(v_minus_inf)));
delta_deg=rad2deg(delta);

[rp, flag] = rpsolver(v_minus_inf, v_plus_inf, flyby.planetId);

% GA altitude
h_GA=rp-astroConstants(25);

%velocities of hiperbolic arcs at pericentre
vp_minus=sqrt((2*astroConstants(15)/rp)+norm(v_minus_inf)^2);
vp_plus=sqrt((2*astroConstants(15)/rp)+norm(v_plus_inf)^2);
delta_vp=abs(vp_minus-vp_plus);

%characterize entry leg hiperbola
e_minus=((rp*+norm(v_minus_inf)^2)/astroConstants(15))+1;
a_minus=-rp/(e_minus-1);


%characterize exit leg hiperbola
e_plus=((rp*+norm(v_plus_inf)^2)/astroConstants(15))+1;
a_plus=-rp/(e_plus-1);

u=[0 0 1];
totalV_vector = v_plus_inf+v_minus_inf;
rp_vector = -rp*(cross(u,totalV_vector))/norm(cross(u,totalV_vector));

%Propagation of orbit minus leg
% Set time span
tspan_minus = linspace( 0, 100*60*60, 5000 );
% Set initial conditions
v_0_minus = sqrt(astroConstants(15)*(2/rp - 1/a_minus));
y0_minus = [rp_vector -v_0_minus.*(totalV_vector./norm(totalV_vector))];
% Perform the integration
[ t_minus, Y_minus ] = ode113( @(t_minus,y_minus) ode_2bp(t_minus,y_minus,astroConstants(15)), tspan_minus, y0_minus, options );


%Propagation of orbit plus leg
% Set time span
tspan_plus = linspace( 0, 100*60*60, 5000 );
% Set initial conditions
v_0_plus = sqrt(astroConstants(15)*(2/rp - 1/a_plus));
y0_plus = [rp_vector v_0_plus.*(totalV_vector./norm(totalV_vector))];
% Perform the integration
[ t_plus, Y_plus ] = ode113( @(t_plus,y_plus) ode_2bp(t_plus,y_plus,astroConstants(15)), tspan_plus, y0_plus, options );


figure()
opts.Units = 'km';
jupiterPlot;
hold on;
plot3( Y_minus(:,1), Y_minus(:,2), Y_minus(:,3),Y_plus(:,1), Y_plus(:,2), Y_plus(:,3));
hold off
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
grid on



%% function find dvmin

function [min, pos1, pos2, pos3] = findMin(dv)
    min = max(max(max(dv)));
    pos1 = 0;
    pos2 = 0;
    pos3 = 0;
    for i = 1:size(dv, 1) 
        for j = 1:size(dv, 2)
            for k = 1:size(dv, 3) 
                if dv(i, j, k) < min
                    min = dv(i, j, k);
                    pos1 = i;
                    pos2 = j;
                    pos3 = k;
                end
            end
        end
    end

    if pos1 == 0 && pos2 == 0 && pos3 == 0
        error("Unable to find minimum")
    end
end

