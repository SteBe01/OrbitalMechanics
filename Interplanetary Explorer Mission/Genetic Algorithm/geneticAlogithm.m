%% Genetic Alogorithm

clear, clc
close all

restoredefaultpath
addpath(genpath("..\\Functions\"))
addpath("..\\Functions_custom\")
addpath("Functions\")

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;
n_iter = 1;
windowChoice = 1;

if windowChoice == 1
    mission.dep_time_lb = [2030 12 09 0 0 0];
    mission.dep_time_ub = [2039 02 25 0 0 0];
    mission.flyby_time_lb = [2042 12 26 0 0 0];
    mission.flyby_time_ub = [2048 09 25 0 0 0];
    mission.arr_time_lb = [2045 06 13 0 0 0];
    mission.arr_time_ub = [2058 01 01 0 0 0];
elseif windowChoice == 2
    mission.dep_time_lb = [2028 01 01 0 0 0];
    mission.dep_time_ub = [2058 01 01 0 0 0];
    mission.flyby_time_lb = [2028 01 01 0 0 0];
    mission.flyby_time_ub = [2058 01 01 0 0 0];
    mission.arr_time_lb = [2028 01 01 0 0 0];
    mission.arr_time_ub = [2058 01 01 0 0 0];
end

lb = [date2mjd2000(mission.dep_time_lb) date2mjd2000(mission.flyby_time_lb) date2mjd2000(mission.arr_time_lb)]';
ub = [date2mjd2000(mission.dep_time_ub) date2mjd2000(mission.flyby_time_ub) date2mjd2000(mission.arr_time_ub)]';

% rng default

options = optimoptions('ga', 'MaxStallGenerations', 15, 'FunctionTolerance', ...
        1e-4, 'MaxGenerations', 1e3, 'NonlinearConstraintAlgorithm', 'penalty',...
        'PopulationSize', 100, 'UseVectorized', false, 'UseParallel',false, InitialPopulationRange=[lb ub]', CrossoverFcn='crossoverlaplace', PlotFcn='gaplotdistance');
if n_iter > 1
    options.Display = "final";
else
    options.Display = "iter";
end

opts = optimset('TolX', eps(1), 'TolFun', eps(1), 'Display', 'off');

for i = 1:n_iter
    [sol,fval,exitflag,output,population,scores] = ga(@(vect) completeInterplanetaryGA(vect(1), vect(2), vect(3), departure.planetId, flyby.planetId, arrival.bodyId), 3, [1 -1 0; 0 1 -1; 1 0 -1], [0 0 0], [], [], lb, ub, [], options);

    [tspan, dv_fmin] = fmincon(@(tspan) completeInterplanetaryGA(tspan(1), tspan(2), tspan(3), departure.planetId, flyby.planetId, arrival.bodyId), [sol(1), sol(2), sol(3)]', [], [], [], [], lb, ub, [], opts);

    [dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(tspan(1), tspan(2), tspan(3), departure.planetId, flyby.planetId, arrival.bodyId);
    
    writematrix([readmatrix("GaResults.csv"); windowChoice, dv_fmin, tspan', dv1, dv2, dv3, rp], "GaResults.csv")
end

disp("Finished!")

