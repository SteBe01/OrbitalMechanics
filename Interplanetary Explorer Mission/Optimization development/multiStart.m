%% Multi Start

% This script was used to create the "Multi Start", later implemented
% as a function and moved to the main script

clear, clc
close all

restoredefaultpath
addpath(genpath("..\\Functions\"))
addpath(genpath("..\\Functions_webeep\"))

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

% rng shuffle

n_elements = 1e4;
parallel = 1;
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

tic

% random vector of players
r1 = zeros(n_elements, 1);
r2 = r1;
r3 = r1;
for i = 1:n_elements
    rand1 = date2mjd2000(mission.dep_time_lb) + (date2mjd2000(mission.dep_time_ub)-date2mjd2000(mission.dep_time_lb)).*rand(1,1);
    rand2 = date2mjd2000(mission.flyby_time_lb) + (date2mjd2000(mission.flyby_time_ub)-date2mjd2000(mission.flyby_time_lb)).*rand(1,1);
    rand3 = date2mjd2000(mission.arr_time_lb) + (date2mjd2000(mission.arr_time_ub)-date2mjd2000(mission.arr_time_lb)).*rand(1,1);

    if rand2 < rand1
        rand2 = rand1 + (date2mjd2000(mission.flyby_time_ub)-rand1).*rand(1,1);
    end
    if rand3 < rand2
        rand3 = rand2 + (date2mjd2000(mission.arr_time_ub)-rand2).*rand(1,1);
    end
    r1(i, :) = rand1;
    r2(i, :) = rand2;
    r3(i, :) = rand3;
end
stpts = [r1, r2, r3];

f = @(vect) completeInterplanetaryMS(vect(1), vect(2), vect(3), departure.planetId, flyby.planetId, arrival.bodyId);
% initial guess
x0(1) = date2mjd2000(mission.dep_time_lb) + (date2mjd2000(mission.dep_time_ub)-date2mjd2000(mission.dep_time_lb)).*rand(1,1);
x0(2) = date2mjd2000(mission.flyby_time_lb) + (date2mjd2000(mission.flyby_time_ub)-date2mjd2000(mission.flyby_time_lb)).*rand(1,1);
x0(3) = date2mjd2000(mission.arr_time_lb) + (date2mjd2000(mission.arr_time_ub)-date2mjd2000(mission.arr_time_lb)).*rand(1,1);
if x0(2) < x0(1)
    x0(2) = x0(1) + (date2mjd2000(mission.flyby_time_ub)-x0(1)).*rand(1,1);
end
if x0(3) < x0(2)
    x0(3) = x0(2) + (date2mjd2000(mission.arr_time_ub)-x0(2)).*rand(1,1);
end

startpts = CustomStartPointSet(stpts);
lb = date2mjd2000(mission.dep_time_lb) * ones(3, 1);
ub = date2mjd2000(mission.arr_time_ub) * ones(3, 1);

opts = optimoptions(@fmincon, 'Algorithm', 'sqp');
newprob = createOptimProblem('fmincon', 'x0', x0, 'lb', lb, 'ub', ub, 'objective', f, 'options', opts);

gs = GlobalSearch;

if parallel
    ms = MultiStart(gs,'UseParallel', true, 'Display', 'iter');
    pool = gcp('nocreate');
    if isempty(pool)
        parpool
    else
        warning("Using existing parpool")
    end
else
    ms = MultiStart(gs,'UseParallel', false, 'Display', 'iter');
end

[xcust, fcust] = run(ms, newprob, startpts);

[dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(xcust(1), xcust(2), xcust(3), departure.planetId, flyby.planetId, arrival.bodyId);
dvTot = dv1 + dv2 + dv3;

time_elapsed = toc;
disp("Time elapsed: " + time_elapsed + " s")
disp("Data - departure Id: " + departure.planetId + ", flyBy Id: " + flyby.planetId + ", arrival Id: " + arrival.bodyId)
disp("Solution: " + dvTot + " km/s, found at " + xcust(1) + " " + xcust(2) + " " + xcust(3))

