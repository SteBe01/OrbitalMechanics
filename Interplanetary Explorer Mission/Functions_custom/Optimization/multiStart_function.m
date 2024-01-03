function [solution] = multiStart_function(mission)

tic

% random vector of players
r1 = zeros(mission.options.n_elements, 1);
r2 = r1;
r3 = r1;
for i = 1:mission.options.n_elements
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

f = @(vect) completeInterplanetaryMS(vect(1), vect(2), vect(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
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

if mission.options.parallel
    ms = MultiStart(gs, 'UseParallel', true, 'Display', 'iter');
    pool = gcp('nocreate');
    if isempty(pool)
        parpool
    else
        warning("Using existing parpool")
    end
else
    ms = MultiStart(gs, 'UseParallel', false, 'Display', 'iter');
end

[xcust, fcust] = run(ms, newprob, startpts);

[dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(xcust(1), xcust(2), xcust(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
dvTot = dv1 + dv2 + dv3;

time_elapsed = toc;
disp("Time elapsed: " + time_elapsed + " s")
disp("Data - departure Id: " + mission.departure_Id + ", flyBy Id: " + mission.flyby_Id + ", arrival Id: " + mission.arrival_Id)
disp("Solution: " + dvTot + " km/s, found at " + xcust(1) + " " + xcust(2) + " " + xcust(3))

solution.dv_fmin = fcust;
solution.tspan = xcust;

end

