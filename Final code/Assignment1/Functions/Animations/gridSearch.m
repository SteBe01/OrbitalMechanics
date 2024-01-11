%% grid search animation

% This script creates the .png files needed to create the animation

clc, clear
close all

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058

restoredefaultpath
addpath(genpath("..\"))

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

% grid search options
mission.options.windowType = 1;
mission.options.window_size = 30;
mission.options.fixedtol = 1e3;                         % in seconds
mission.options.fmincon_choice = 3;                     % 0 for no fmincon
mission.options.animation = 1;

[results] = gridSearch_function_animation(mission);

