%% contour plots (with or without dv min position and with or without time window): departure - flyby and flyby - arrival

% This script prints the contour plots and the synodic periods we used to
% determine the final time window for the optimizations (the two scripts
% are independent from eachother, move to line 187 to run the synodic
% period).
% Choose the parameters for the contour plots between lines 13 and 15.

clc, clear
close all

% ------------------------------------------------
toggle_timeWindow = 1;
toggle_finalSolution = 1;
toggle_titles = 1;
% ------------------------------------------------

restoredefaultpath
addpath(genpath("."))

% our solution
load("dvmin.mat")

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 500);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 500);

% contour plots - departure, fly by
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
        [flyby.kep, ~] = uplanet(arr_time_vect(j), flyby.planetId);

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

% dv = dv_1 + dv_2;
dv = dv_1; % without flyby dv

plot1 = figure();
contour(dep_time_vect, arr_time_vect, dv', 2:0.15:8, HandleVisibility="off")
if toggle_titles
    title("Contour plot: departure - flyby")
end
c = colorbar;
c.Label.String = '[km/s]';
grid on, hold on
xlabel("Departure [mjd2000]")
ylabel("Flyby [mjd2000]")

if toggle_finalSolution
    [pos1, pos2] = find(dv == min(min(dv)));
    plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1.7, MarkerSize=8)
    plot(xcust(1), xcust(2), 'or', LineWidth=1.7, MarkerSize=6)

    if ~toggle_timeWindow
        legend("Contour plot min", "Mission min", Location="southeast")
    end
end
% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")

% contour plots - fly by, arrival
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

% dv = dv_1 + dv_2;
dv = dv_2; % without flyby dv

plot2 = figure();
contour(dep_time_vect, arr_time_vect, dv', 4:0.5:20, HandleVisibility="off")
if toggle_titles
    title("Contour plot: flyby - arrival")
end
c = colorbar;
c.Label.String = '[km/s]';
grid on, hold on, axis equal
xlabel("Flyby [mjd2000]")
ylabel("Arrival [mjd2000]")

if toggle_finalSolution
    [pos1, pos2] = find(dv == min(min(dv)));
    plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1.7, MarkerSize=8)
    plot(xcust(2), xcust(3), 'or', LineWidth=1.7, MarkerSize=6)

    if ~toggle_timeWindow
        legend("Contour plot min", "Mission min", Location="southeast")
    end
end
% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")

if toggle_timeWindow
    mission.dep_time_lb = date2mjd2000([2030 12 09 0 0 0]);
    mission.dep_time_ub = date2mjd2000([2039 02 25 0 0 0]);
    mission.flyby_time_lb = date2mjd2000([2042 12 26 0 0 0]);
    mission.flyby_time_ub = date2mjd2000([2048 09 25 0 0 0]);
    mission.arr_time_lb = date2mjd2000([2045 06 13 0 0 0]);
    mission.arr_time_ub = date2mjd2000([2058 01 01 0 0 0]);
    
    figure(plot1)
    plot([mission.dep_time_lb mission.dep_time_ub], [mission.flyby_time_lb mission.flyby_time_lb], 'r', LineWidth=2, HandleVisibility='off');
    plot([mission.dep_time_lb mission.dep_time_ub], [mission.flyby_time_ub mission.flyby_time_ub], 'r', LineWidth=2, HandleVisibility='off');
    plot([mission.dep_time_lb mission.dep_time_lb], [mission.flyby_time_lb mission.flyby_time_ub], 'r', LineWidth=2, HandleVisibility='off');
    if toggle_finalSolution
        plot([mission.dep_time_ub mission.dep_time_ub], [mission.flyby_time_lb mission.flyby_time_ub], 'r', LineWidth=2);
        legend("Contour plot min", "Mission min", "Time window", Location="southeast")
    else
        plot([mission.dep_time_ub mission.dep_time_ub], [mission.flyby_time_lb mission.flyby_time_ub], 'r', LineWidth=2, HandleVisibility='off');
    end
    
    figure(plot2)
    plot([mission.flyby_time_lb mission.flyby_time_ub], [mission.arr_time_lb mission.arr_time_lb], 'r', LineWidth=2, HandleVisibility='off');
    plot([mission.flyby_time_lb mission.flyby_time_ub], [mission.arr_time_ub mission.arr_time_ub], 'r', LineWidth=2, HandleVisibility='off');
    plot([mission.flyby_time_lb mission.flyby_time_lb], [mission.arr_time_lb mission.arr_time_ub], 'r', LineWidth=2, HandleVisibility='off');
    if toggle_finalSolution
        plot([mission.flyby_time_ub mission.flyby_time_ub], [mission.arr_time_lb mission.arr_time_ub], 'r', LineWidth=2);
        legend("Contour plot min", "Mission min", "Time window", Location="southeast")
    else
        plot([mission.flyby_time_ub mission.flyby_time_ub], [mission.arr_time_lb mission.arr_time_ub], 'r', LineWidth=2, HandleVisibility='off');
    end
end

fontsize(plot1, 15, "points")
fontsize(plot2, 15, "points")


%% Synodic Periods

clc, clear

restoredefaultpath
addpath(genpath("."))

date_departure = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(date_departure);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.planetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);

departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [s]
flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [s]
arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [s]

%Synodic Period between departure planet and flyby
T_syn_dep_flyby=(departure.T_orb*flyby.T_orb)/(abs(departure.T_orb-flyby.T_orb))/(60*60*24*astroConstants(32)); %[years]

%Synodic Period between flyby planet and arrival
T_syn_flyby_arr=(arrival.T_orb*flyby.T_orb)/(abs(arrival.T_orb-flyby.T_orb))/(60*60*24*astroConstants(32)); %[years]

%Synodic Period between departure planet and arrival
T_syn_dep_arr=(arrival.T_orb*departure.T_orb)/(abs(arrival.T_orb-departure.T_orb))/(60*60*24*astroConstants(32)); %[years]

fprintf("Synodic period between Departure and Flyby: \t" + T_syn_dep_flyby + " years\n");
fprintf("Synodic period between Flyby and Arrival: \t\t" + T_syn_flyby_arr + " years\n");

