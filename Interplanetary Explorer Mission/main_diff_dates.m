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

%% porkchop plots - departure, fly by

clc, clear
close all

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

departure.planetId = 6;
flyby.plnetId = 5;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);

% dep_time_vect = linspace(date2mjd2000(dep_time), 1.45e4, 100);
% arr_time_vect = linspace(1.47e4, 1.78e4, 100);

orbitType = 0;
dv_1 = zeros(length(dep_time_vect), length(arr_time_vect));
dv_2 = dv_1;
tof_vect = dv_1;
for i = 1:length(dep_time_vect)
    for j = 1:length(arr_time_vect)
        tof = (arr_time_vect(j) - dep_time_vect(i)) * 24 * 60 * 60; % seconds

        if tof <= 1e5
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end
        tof_vect(i, j) = tof;

        [departure.kep, ksun] = uplanet(dep_time_vect(i), departure.planetId);
        [flyby.kep, ~] = uplanet(arr_time_vect(j), flyby.plnetId);

        [departure.r0, departure.v0] = kep2car([departure.kep ksun]);
        [flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(departure.r0, flyby.r0, tof, ksun, orbitType, 0);
        if A < 0
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end

        dv_1(i, j) = norm(VI - departure.v0);
        dv_2(i, j) = norm(flyby.v0 - VF);
    end
end

dv = dv_1 + dv_2;

contour(dep_time_vect, arr_time_vect, dv', 2:0.2:8)
colorbar, grid on, hold on
xlabel("Departure")
ylabel("Arrival")

[pos1, pos2] = find(dv == min(min(dv)));
plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1)

% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")



%% porkchop plots - fly by, arrival

clc, clear
close all

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

flyby.planetId = 5;
arrival.bodyId = 79;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);

% dep_time_vect = linspace(1.35e4, 1.48e4, 100);
% arr_time_vect = linspace(1.535e4, 1.565e4, 100);

orbitType = 0;
dv_1 = zeros(length(dep_time_vect), length(arr_time_vect));
dv_2 = dv_1;
tof_vect = dv_1;
for i = 1:length(dep_time_vect)
    for j = 1:length(arr_time_vect)
        tof = (arr_time_vect(j) - dep_time_vect(i)) * 24 * 60 * 60; % seconds

        if tof <= 1e5
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end
        tof_vect(i, j) = tof;

        [flyby.kep, ksun] = uplanet(dep_time_vect(i), flyby.planetId);
        [arrival.kep, ~] = ephNEO(arr_time_vect(j), arrival.bodyId);

        [flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);
        [arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(flyby.r0, arrival.r0, tof, ksun, orbitType, 0);
        if A < 0
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end

        dv_1(i, j) = norm(VI - flyby.v0);
        dv_2(i, j) = norm(arrival.v0 - VF);
    end
end

dv = dv_1 + dv_2;

contour(dep_time_vect, arr_time_vect, dv', 12:0.5:25) % 12:35
colorbar, grid on, hold on, axis equal
xlabel("Departure")
ylabel("Arrival")

[pos1, pos2] = find(dv == min(min(dv)));
plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1)

% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")


%% Grid Search for departure-flyby-arrival 
clc


%--------Define well this values with Tsyn and ToF!!!!!!
mission_dep_time = date_departure;
mission_flyby_time=date_flyby;
mission_arr_time = date_arrival;
%-----------------------------------------------------

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;
fixedtol = 1e3;
window_size = 20;

time1 = date2mjd2000(mission_dep_time);
time2=date2mjd2000(mission_flyby_time);
time3 = date2mjd2000(mission_arr_time);


departure.time_vect = linspace(time1, time2, window_size);
flyby.time_vect = linspace(time2-5, time2+5, window_size); %time required for flyby basically
arrival.time_vect = linspace(time2, time3, window_size);

orbitType = 0;
dv_1 = zeros(length(departure.time_vect), length(flyby.time_vect), length(arrival.time_vect));
dv_2 = dv_1;
dv_3 = dv_1;
dv_flyby_tot = dv_1;
tof_vect_1 = dv_1;
tof_vect_2 = dv_1;
num_iter = 0;

tol = (departure.time_vect(end)-departure.time_vect(1)) * 24 * 60 * 60;

while fixedtol < tol
    tic
    num_iter = num_iter + 1;
    iteration.number{num_iter} = num_iter;

    for i = 1:length(departure.time_vect)
        for j = 1:length(flyby.time_vect)
            for k = 1:length(arrival.time_vect)
    
                tof_1 = (flyby.time_vect(j) - departure.time_vect(i)) * 24 * 60 * 60; % seconds
                tof_2 = (arrival.time_vect(k) - flyby.time_vect(j)) * 24 * 60 * 60; % seconds
        
                if (tof_1 <= 1e5) && (tof_2 <= 1e5)
                    dv_1(i, j, k) = NaN;
                    dv_2(i, j, k) = NaN;
                    dv_3(i, j, k) = NaN;
                    continue
                end
                tof_vect_1(i, j, k) = tof_1;
                tof_vect_2(i, j, k) = tof_2;
        
                [departure.kep, ksun] = uplanet(departure.time_vect(i), departure.planetId);
                [flyby.kep, flyby.mu] = uplanet(flyby.time_vect(j), flyby.planetId);
                [arrival.kep, ~] = ephNEO(arrival.time_vect(k), arrival.bodyId);
        
                [departure.r0, departure.v0] = kep2car([departure.kep ksun]);
                [flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);
                [arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);
        
                [A_1, P_1, E_1, ERROR_1, VI_1, VF_1, TPAR_1, THETA_1] = lambertMR(departure.r0, flyby.r0, tof_1, ksun, orbitType, 0);
                if A_1 < 0
                    dv_1(i, j, k) = NaN;
                    dv_2(i, j, k) = NaN;
                    dv_3(i, j, k) = NaN;
                    continue
                end
    
                [A_2, P_2, E_2, ERROR_2, VI_2, VF_2, TPAR_2, THETA_2] = lambertMR(flyby.r0, arrival.r0, tof_2, ksun, orbitType, 0);
                if A_2 < 0
                    dv_1(i, j, k) = NaN;
                    dv_2(i, j, k) = NaN;
                    dv_3(i, j, k) = NaN;
                    continue
                end
        
                v_inf_minus = VI_2 + flyby.v0;
                v_inf_plus = VF_1 + flyby.v0;
                dv_1(i, j, k) = norm(VI_1 - departure.v0);
                [rp, flag] = rpsolver(v_inf_minus, v_inf_plus, flyby.planetId);
                if flag == 0
                    dv_1(i, j, k) = NaN;
                    dv_2(i, j, k) = NaN;
                    dv_3(i, j, k) = NaN;
                    continue
                end
                dv_2(i, j, k) = abs(sqrt((2*astroConstants(flyby.planetId + 10)/rp)+norm(v_inf_plus)^2)-sqrt((2*astroConstants(flyby.planetId + 10)/rp)+norm(v_inf_minus)^2));
                dv_3(i, j, k) = norm(arrival.v0 - VF_2);
                dv_flyby_tot(i, j, k) = norm(VF_1 - VF_2);
            end
        end
    end

    dv = dv_1 + dv_2 + dv_3;
    [dVmin, pos1, pos2, pos3] = findMin(dv);

    if pos1 == 1
        pos1_bnd_1 = 1;
        pos1_bnd_2 = 1 + 1;
    elseif pos1 == length(departure.time_vect)
        pos1_bnd_2 = length(departure.time_vect);
        pos1_bnd_1 = length(departure.time_vect) - 1;
    else
        pos1_bnd_1 = pos1 - 1;
        pos1_bnd_2 = pos1 + 1;
    end

    if pos2 == 1
        pos2_bnd_1 = 1;
        pos2_bnd_2 = 1 + 1;
    elseif pos2 == length(flyby.time_vect)
        pos2_bnd_2 = length(flyby.time_vect);
        pos2_bnd_1 = length(flyby.time_vect) - 1;
    else
        pos2_bnd_1 = pos2 - 1;
        pos2_bnd_2 = pos2 + 1;
    end

    if pos3 == 1
        pos3_bnd_1 = 1;
        pos3_bnd_2 = 1 + 1;
    elseif pos3 == length(arrival.time_vect)
        pos3_bnd_2 = length(arrival.time_vect);
        pos3_bnd_1 = length(arrival.time_vect) - 1;
    else
        pos3_bnd_1 = pos3 - 1;
        pos3_bnd_2 = pos3 + 1;
    end
    tol = (departure.time_vect(pos1_bnd_2) - departure.time_vect(pos1_bnd_1)) * 24 * 60 * 60;

    if tol > fixedtol
        departure.time_vect = linspace(departure.time_vect(pos1_bnd_1), departure.time_vect(pos1_bnd_2), window_size);
        flyby.time_vect = linspace(flyby.time_vect(pos2_bnd_1), flyby.time_vect(pos2_bnd_2), window_size);
        arrival.time_vect = linspace(arrival.time_vect(pos3_bnd_1), arrival.time_vect(pos3_bnd_2), window_size);
    end

    iteration.time{num_iter} = toc;
    iteration.dv1{num_iter} = dv_1;
    iteration.dv2{num_iter} = dv_2;
    iteration.dv3{num_iter} = dv_3;
    iteration.dv{num_iter} = dv;
    iteration.dv_min{num_iter} = dVmin;
    iteration.dv_2_min{num_iter} = min(min(min(dv_2)));
    disp ("Iteration " + num_iter + " done in " + toc + " s, dv = " + dVmin + " km/s")
end

%-PLOT-----------

% Define Actual Parameters

departure.Date = mjd20002date(departure.time_vect(pos1));
flyby.Date = mjd20002date(flyby.time_vect(pos2));
ToF_dep_flyby=(flyby.time_vect(pos2)-departure.time_vect(pos1))*24*60*60;

arrival.Date = mjd20002date(arrival.time_vect(pos3));
ToF_flyby_arr=(arrival.time_vect(pos3)-flyby.time_vect(pos2))*24*60*60;

[departure.kep_actual, ksun_actual] = uplanet(departure.time_vect(pos1), departure.planetId);
[flyby.kep_actual, ~] = uplanet(flyby.time_vect(pos2), flyby.planetId);
[arrival.kep_actual, ~] = ephNEO(arrival.time_vect(pos3), arrival.bodyId);

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% orbit propagation - departure
[departure.r0_actual, departure.v0_actual] = kep2car([departure.kep_actual, ksun_actual]);
departure.y0_actual = [departure.r0_actual departure.v0_actual];

departure.T_orb_actual = 2*pi*sqrt( departure.kep_actual(1)^3/ksun_actual ); % Orbital period [1/s]
departure.tspan= linspace( 0, departure.T_orb_actual, 200 );
[ ~, departure.Y_actual ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), departure.tspan, departure.y0_actual, options );

% orbit propagation - flyby
[flyby.r0_actual, flyby.v0_actual] = kep2car([flyby.kep_actual, ksun_actual]);
flyby.y0_actual = [flyby.r0_actual flyby.v0_actual];

flyby.T_orb_actual = 2*pi*sqrt( flyby.kep_actual(1)^3/ksun_actual ); % Orbital period [1/s]
flyby.tspan= linspace( 0, flyby.T_orb_actual, 200 );
[ ~, flyby.Y_actual ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), flyby.tspan, flyby.y0_actual, options );

% orbit propagation - arrival
[arrival.r0_actual, arrival.v0_actual] = kep2car([arrival.kep_actual, ksun_actual]);
arrival.y0_actual = [arrival.r0_actual arrival.v0_actual];

arrival.T_orb_actual = 2*pi*sqrt( arrival.kep_actual(1)^3/ksun_actual ); % Orbital period [1/s]
arrival.tspan= linspace( 0, arrival.T_orb_actual, 200 );
[ ~, arrival.Y_actual ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), arrival.tspan, arrival.y0_actual, options );

%Propagation first transfer ARC
[A_1_actual, P_1_actual, E_1_actual, ERROR_1_actual, VI_1_actual, VF_1_actual, TPAR_1_actual, THETA_1_actual] = lambertMR(departure.r0_actual, flyby.r0_actual, ToF_dep_flyby, ksun_actual, orbitType, 0);
y0_1 = [ departure.r0_actual VI_1_actual ];
% Set time span
T_1 = 2*pi*sqrt( A_1_actual^3/ksun_actual ); % Orbital period [s]
tspan_1 = linspace( 0,ToF_dep_flyby, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% Perform the integration
[   t, Y_1 ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), tspan_1, y0_1, options );

%Propagation second transfer ARC
[A_2_actual, P_2_actual, E_2_actual, ERROR_2_actual, VI_2_actual, VF_2_actual, TPAR_2_actual, THETA_2_actual] = lambertMR(flyby.r0_actual, arrival.r0_actual, ToF_flyby_arr, ksun_actual, orbitType, 0);
y0_2 = [ flyby.r0_actual VI_2_actual ];
% Set time span
T_2 = 2*pi*sqrt( A_2_actual^3/ksun_actual ); % Orbital period [s]
tspan_2 = linspace( 0,ToF_flyby_arr, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% Perform the integration
[   t, Y_2 ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), tspan_2, y0_2, options );

% Plot the results
figure()
plot3( departure.Y_actual(:,1), departure.Y_actual(:,2), departure.Y_actual(:,3), '-','color', 'b' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('orbits');
axis equal;
grid on;
hold on
plot3( flyby.Y_actual(:,1),flyby.Y_actual(:,2), flyby.Y_actual(:,3), '-' ,'color', 'r')
plot3( arrival.Y_actual(:,1), arrival.Y_actual(:,2), arrival.Y_actual(:,3), '-','color', 'g')
plot3( Y_1(:,1), Y_1(:,2), Y_1(:,3), '--','color', 'm' )
plot3( Y_2(:,1), Y_2(:,2), Y_2(:,3), '--','color', 'c' )
plot3(departure.r0_actual(1),departure.r0_actual(2),departure.r0_actual(3),'o','Color','b','MarkerFaceColor','b')
plot3(flyby.r0_actual(1),flyby.r0_actual(2),flyby.r0_actual(3),'o','Color','r','MarkerFaceColor','r')
plot3(arrival.r0_actual(1),arrival.r0_actual(2),arrival.r0_actual(3),'o','Color','g','MarkerFaceColor','g')
plot3(0,0,0,'o','Color','y','MarkerFaceColor','y')
legend('Saturn Orbit','Jupiter Orbit','Asteroid N.79 Orbit','Transfer Arc 1','Transfer Arc 2', ...
    'Saturn','Jupiter','Asteroid N.79','Sun')
hold off


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
