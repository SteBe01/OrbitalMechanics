%% Multi Start

clear, clc
close all

restoredefaultpath
addpath(genpath("..\\Functions\"))
addpath("..\\Functions_custom\")
addpath("Functions\")

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

n_elements = 1e4;
parallel = 1;

mission.dep_time = [2028 01 01 0 0 0];
mission.arr_time = [2058 01 01 0 0 0];

tic

% random vector of players
r1 = zeros(n_elements, 1);
r2 = r1;
r3 = r1;
for i = 1:n_elements
    rand1 = date2mjd2000(mission.dep_time) + (date2mjd2000(mission.arr_time)-date2mjd2000(mission.dep_time)).*rand(1,1);
    rand2 = rand1 + (date2mjd2000(mission.arr_time)-rand1).*rand(1,1);
    rand3 = rand2 + (date2mjd2000(mission.arr_time)-rand2).*rand(1,1);
    r1(i, :) = rand1;
    r2(i, :) = rand2;
    r3(i, :) = rand3;
end
stpts = [r1, r2, r3];

f = @(vect) completeInterplanetaryMS(vect(1), vect(2), vect(3), departure.planetId, flyby.planetId, arrival.bodyId);
% initial guess
x0(1) = date2mjd2000(mission.dep_time) + (date2mjd2000(mission.arr_time)-date2mjd2000(mission.dep_time)).*rand(1,1);
x0(2) = x0(1) + (date2mjd2000(mission.arr_time)-x0(1)).*rand(1,1);
x0(3) = x0(2) + (date2mjd2000(mission.arr_time)-x0(2)).*rand(1,1);

startpts = CustomStartPointSet(stpts);
lb = date2mjd2000(mission.dep_time) * ones(3, 1);
ub = date2mjd2000(mission.arr_time) * ones(3, 1);

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

