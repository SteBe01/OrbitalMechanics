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

mission.dep_time = [2028 01 01 0 0 0];
mission.arr_time = [2058 01 01 0 0 0];

f = @(vect) completeInterplanetaryMS(vect(1), vect(2), vect(3), departure.planetId, flyby.planetId, arrival.bodyId);

r1 = date2mjd2000(mission.dep_time) + (date2mjd2000(mission.arr_time)-date2mjd2000(mission.dep_time)).*rand(n_elements,1);
r2 = date2mjd2000(mission.dep_time) + (date2mjd2000(mission.arr_time)-date2mjd2000(mission.dep_time)).*rand(n_elements,1);
r3 = date2mjd2000(mission.dep_time) + (date2mjd2000(mission.arr_time)-date2mjd2000(mission.dep_time)).*rand(n_elements,1);
x0 = date2mjd2000(mission.dep_time) + (date2mjd2000(mission.arr_time)-date2mjd2000(mission.dep_time)).*rand(3,1);

stpts = [r1, r2, r3];
startpts = CustomStartPointSet(stpts);
lb = date2mjd2000(mission.dep_time);
ub = date2mjd2000(mission.arr_time);

opts = optimoptions(@fminunc, 'Algorithm', 'quasi-newton');
newprob = createOptimProblem('fminunc', 'x0', x0, 'lb', lb, 'ub', ub, 'objective', f, 'options', opts);

gs = GlobalSearch("Display","iter");
ms = MultiStart(gs);

[xcust, fcust] = run(ms, newprob, startpts);

[dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(xcust(1), xcust(2), xcust(3), departure.planetId, flyby.planetId, arrival.bodyId);
dvTot = dv1 + dv2 + dv3;

