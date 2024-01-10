%% data

% This program is the main script of our Interplanetary Mission.
% By specifing the simChoice (line 14) it can use three different
% optimization techniques (1: Grid Search, 2: Genetic Algorithm and
% 3: Multi Start).
% It runs automatically and prints all the important data relative to the
% mission. Multiple plots will be opened and an animation will play (only
% with "Grid Search").

clc, clear
close all

simChoice = 1;

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058

restoredefaultpath
addpath(genpath("Functions\"))
addpath(genpath("Functions_webeep\"))

mission.departure_Id = 6;
mission.flyby_Id = 5;
mission.arrival_Id = 79;

mission.dep_time_lb = [2030 12 09 0 0 0];
mission.dep_time_ub = [2039 02 25 0 0 0];
mission.flyby_time_lb = [2042 12 26 0 0 0];
mission.flyby_time_ub = [2048 09 25 0 0 0];
mission.arr_time_lb = [2045 06 13 0 0 0];
mission.arr_time_ub = [2058 01 01 0 0 0];

mission.dep_time = [2028 01 01 0 0 0];
mission.arr_time = [2058 01 01 0 0 0];

startInfoDisplay(mission)
windowDisplay(mission)

if simChoice == 1
    disp('"Grid Search" optimization started...')

    % grid search options
    mission.options.windowType = 1;
    mission.options.window_size = 30;
    mission.options.fixedtol = 1e3;                         % in seconds
    mission.options.fmincon_choice = 3;                     % 0 for no fmincon
    mission.options.animation = 1;

    [results] = gridSearch_function(mission);
elseif simChoice == 2
    disp('"Genetic Algorithm" optimization started...')

    % ga options
    mission.options.n_iter = 5;

    [results] = ga_function(mission);
elseif simChoice == 3
    disp('"Multi Start" optimization started...')

    % multi start options
    mission.options.n_elements = 1e4;
    mission.options.parallel = 1;

    [results] = multiStart_function(mission);
else
    error("Invalid simulation choice! Try 1, 2 or 3")
end

missionPlot(results.tspan(1), results.tspan(2), results.tspan(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
[data, ~, ~] = flybyPlot(results.tspan(1), results.tspan(2), results.tspan(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id, 1e7);

[results.dv1, results.dv2, results.dv3, results.rp, results.exitValue] = completeInterplanetary(results.tspan(1), results.tspan(2), results.tspan(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
endInfoDisplay(data, results)

