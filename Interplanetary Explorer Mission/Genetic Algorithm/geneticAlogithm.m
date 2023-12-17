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

time_instant_dep = [2028 01 01 0 0 0];
time_instant_arr = [2058 01 01 0 0 0];
lb = ones(3,1) .* date2mjd2000(time_instant_dep);
ub = ones(3,1) .* date2mjd2000(time_instant_arr);

options = optimoptions('ga', 'MaxStallGenerations', 15, 'FunctionTolerance', ...
        1e-4, 'MaxGenerations', 100, 'NonlinearConstraintAlgorithm', 'penalty',...
        'PopulationSize', 100, 'Display', 'iter', 'UseVectorized', false, 'UseParallel',false, 'MutationFcn','mutationadaptfeasible');

[sol,fval,exitflag,output,population,scores] = ga(@(vect) completeInterplanetaryGA(vect(1), vect(2), vect(3), departure.planetId, flyby.planetId, arrival.bodyId), 3, [1 -1 0; 0 1 -1; 1 0 -1], [0 0 0], [], [], lb, ub, [], options);

[dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(sol(1), sol(2), sol(3), departure.planetId, flyby.planetId, arrival.bodyId);

