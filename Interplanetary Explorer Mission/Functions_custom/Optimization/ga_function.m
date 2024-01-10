function [solution] = ga_function(mission)

% Core function for genetic algorithm optimization
%
% Usage
% [solution] = ga_function(mission)
%
% Input arguments:
% ----------------------------------------------------------------
% mission       [-]       mission data          [struct]
%
% Output arguments:
% -----------------------------------------------------------------
% solution      [-]       mission solution      [struct]
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

lb = [date2mjd2000(mission.dep_time_lb) date2mjd2000(mission.flyby_time_lb) date2mjd2000(mission.arr_time_lb)]';
ub = [date2mjd2000(mission.dep_time_ub) date2mjd2000(mission.flyby_time_ub) date2mjd2000(mission.arr_time_ub)]';

% rng default

options = optimoptions('ga', 'MaxStallGenerations', 15, 'FunctionTolerance', ...
        1e-4, 'MaxGenerations', 1e3, 'NonlinearConstraintAlgorithm', 'penalty',...
        'PopulationSize', 100, 'UseVectorized', false, 'UseParallel',false, InitialPopulationRange=[lb ub]', CrossoverFcn='crossoverlaplace', PlotFcn='gaplotdistance');
if mission.options.n_iter > 1
    options.Display = "final";
else
    options.Display = "iter";
end

opts = optimset('TolX', eps(1), 'TolFun', eps(1), 'Display', 'off');

sol_temp = zeros(mission.options.n_iter, 4);
for i = 1:mission.options.n_iter
    [sol,fval,exitflag,solutionput,population,scores] = ga(@(vect) completeInterplanetaryGA(vect(1), vect(2), vect(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id), 3, [1 -1 0; 0 1 -1; 1 0 -1], [0 0 0], [], [], lb, ub, [], options);

    [tspan, dv_fmin] = fmincon(@(tspan) completeInterplanetaryGA(tspan(1), tspan(2), tspan(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id), [sol(1), sol(2), sol(3)]', [], [], [], [], lb, ub, [], opts);

    [dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(tspan(1), tspan(2), tspan(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
    
    sol_temp(i, 1) = dv_fmin;
    sol_temp(i, 2:4) = tspan;
    % writematrix([readmatrix("GaResults.csv"); windowChoice, dv_fmin, tspan', dv1, dv2, dv3, rp], "GaResults.csv")
end

solution.dvMin = dv_fmin;
solution.tspan = tspan;

end

