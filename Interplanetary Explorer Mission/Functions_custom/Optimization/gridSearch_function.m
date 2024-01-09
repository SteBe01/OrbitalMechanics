function [solutions] = gridSearch_function(mission)

% Core function for grid search optimization
%
% Usage
% [solutions] = gridSearch_function(mission)
%
% Input arguments:
% ----------------------------------------------------------------
% mission       [-]       mission data          [struct]
%
% -----------------------------------------------------------------
% Output arguments:
% 
% solution      [-]       mission solution      [struct]

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

[departureTime, flybyTime, arrivalTime] = loadWindows(mission.options.windowType);

if mission.options.animation
    porkchop_start(dep_time, arr_time, mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
end

solutions_dvMin = zeros(length(departureTime), 1);
solutions_tspan = zeros(length(departureTime), 3);

for totWindows = 1:length(departureTime)
    departure.time_vect = linspace(departureTime{totWindows}.lb, departureTime{totWindows}.ub, mission.options.window_size);
    flyby.time_vect = linspace(flybyTime{totWindows}.lb, flybyTime{totWindows}.ub, mission.options.window_size);
    arrival.time_vect = linspace(arrivalTime{totWindows}.lb, arrivalTime{totWindows}.ub, mission.options.window_size);
    
    dv_1 = NaN* ones(length(departure.time_vect), length(flyby.time_vect), length(arrival.time_vect));
    dv_2 = dv_1;
    dv_3 = dv_1;
    rp = dv_1;
    num_iter = 0;
    
    tol = (departure.time_vect(end)-departure.time_vect(1)) * 24 * 60 * 60;
    
    while mission.options.fixedtol < tol
        tic
        num_iter = num_iter + 1;
        
        reverseStr = '';
        for i = 1:length(departure.time_vect)
    
            msg = sprintf('Window %i processed %d percent', totWindows, ceil((i / length(departure.time_vect)) * 100));
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
            for j = i:length(flyby.time_vect)
                for k = j:length(arrival.time_vect)
        
                    tof_1 = (flyby.time_vect(j) - departure.time_vect(i)) * 24 * 60 * 60; % seconds
                    tof_2 = (arrival.time_vect(k) - flyby.time_vect(j)) * 24 * 60 * 60; % seconds
            
                    if (tof_1 <= 1e5) && (tof_2 <= 1e5)
                        continue
                    end
            
                    [dv_1(i, j, k), dv_2(i, j, k), dv_3(i, j, k), rp_temp, exitValue] = completeInterplanetary(departure.time_vect(i), flyby.time_vect(j), arrival.time_vect(k), mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
                    if exitValue
                        continue
                    end
                    rp(i, j, k) = rp_temp;
                end
            end
        end
    
        fprintf(": ");
    
        dv = dv_1 + dv_2 + dv_3;
        [dVmin, pos1, pos2, pos3] = findMin3(dv);
    
        if pos1 == 1
            pos1_bnd_1 = 1;
            pos1_bnd_2 = 1 + 1;
        elseif pos1 == length(departure.time_vect)
            pos1_bnd_2 = length(departure.time_vect);
            pos1_bnd_1 = length(departure.time_vect) - 1;
        else
            pos1_bnd_1 = pos1 - 1;
            pos1_bnd_2 = pos1 + 1;
        end
    
        if pos2 == 1
            pos2_bnd_1 = 1;
            pos2_bnd_2 = 1 + 1;
        elseif pos2 == length(flyby.time_vect)
            pos2_bnd_2 = length(flyby.time_vect);
            pos2_bnd_1 = length(flyby.time_vect) - 1;
        else
            pos2_bnd_1 = pos2 - 1;
            pos2_bnd_2 = pos2 + 1;
        end
    
        if pos3 == 1
            pos3_bnd_1 = 1;
            pos3_bnd_2 = 1 + 1;
        elseif pos3 == length(arrival.time_vect)
            pos3_bnd_2 = length(arrival.time_vect);
            pos3_bnd_1 = length(arrival.time_vect) - 1;
        else
            pos3_bnd_1 = pos3 - 1;
            pos3_bnd_2 = pos3 + 1;
        end
        tol = (departure.time_vect(pos1_bnd_2) - departure.time_vect(pos1_bnd_1)) * 24 * 60 * 60;
    
        if tol > mission.options.fixedtol
            departure.time_vect = linspace(departure.time_vect(pos1_bnd_1), departure.time_vect(pos1_bnd_2), mission.options.window_size);
            flyby.time_vect = linspace(flyby.time_vect(pos2_bnd_1), flyby.time_vect(pos2_bnd_2), mission.options.window_size);
            arrival.time_vect = linspace(arrival.time_vect(pos3_bnd_1), arrival.time_vect(pos3_bnd_2), mission.options.window_size);
        end
    
        if mission.options.fmincon_choice == num_iter
            break
        else
            disp ("iteration " + num_iter + " done in " + toc + " s, dv = " + dVmin + " km/s")
        end
    end
    
    if mission.options.fmincon_choice ~= 0
        opts = optimset('TolX', eps(1), 'TolFun', eps(1), 'Display', 'off');
        lb = date2mjd2000(dep_time) * ones(3,1);
        ub = date2mjd2000(arr_time) * ones(3,1);
        [tspan, dv_fmin] = fmincon(@(tspan) completeInterplanetaryGS(tspan(1), tspan(2), tspan(3), mission.departure_Id, mission.flyby_Id, mission.arrival_Id), [departure.time_vect(pos1) flyby.time_vect(pos2) arrival.time_vect(pos3)]', [], [], [], [], lb, ub, [], opts);
        
        disp("iteration " + num_iter + " done in " + toc + " s, dv = " + dv_fmin + " km/s")
    end

    if totWindows > 1
        delete(a1)
        delete(a2)
        delete(a3)
        delete(a4)
        delete(a5)
        delete(b1)
        delete(b2)
        delete(b3)
        delete(b4)
        delete(b5)
    end

    subplot(1, 2, 1)
    a1 = plot([departureTime{totWindows}.lb departureTime{totWindows}.ub], [flybyTime{totWindows}.lb flybyTime{totWindows}.lb], 'r');
    hold on
    a2 = plot([departureTime{totWindows}.lb departureTime{totWindows}.ub], [flybyTime{totWindows}.ub flybyTime{totWindows}.ub], 'r');
    a3 = plot([departureTime{totWindows}.lb departureTime{totWindows}.lb], [flybyTime{totWindows}.lb flybyTime{totWindows}.ub], 'r');
    a4 = plot([departureTime{totWindows}.ub departureTime{totWindows}.ub], [flybyTime{totWindows}.lb flybyTime{totWindows}.ub], 'r');
    a5 = plot(tspan(1), tspan(2), 'xr', LineWidth=1);
    subplot(1, 2, 2)
    b1 = plot([flybyTime{totWindows}.lb flybyTime{totWindows}.ub], [arrivalTime{totWindows}.lb arrivalTime{totWindows}.lb], 'r');
    hold on
    b2 = plot([flybyTime{totWindows}.lb flybyTime{totWindows}.ub], [arrivalTime{totWindows}.ub arrivalTime{totWindows}.ub], 'r');
    b3 = plot([flybyTime{totWindows}.lb flybyTime{totWindows}.lb], [arrivalTime{totWindows}.lb arrivalTime{totWindows}.ub], 'r');
    b4 = plot([flybyTime{totWindows}.ub flybyTime{totWindows}.ub], [arrivalTime{totWindows}.lb arrivalTime{totWindows}.ub], 'r');
    b5 = plot(tspan(2), tspan(3), 'xr', LineWidth=1);
    sgtitle("Delta velocity = " + dv_fmin + " km/s")
    drawnow

    solutions_dvMin(totWindows) = dv_fmin;
    solutions_tspan(totWindows, :) = tspan;
end

[solutions.dvMin, min_dv_index] = min(solutions_dvMin);
solutions.tspan = solutions_tspan(min_dv_index, :);


function [min, pos1, pos2, pos3] = findMin3(dv)
    min = max(max(max(dv)));
    pos1 = 0;
    pos2 = 0;
    pos3 = 0;
    for ii = 1:size(dv, 1) 
        for jj = 1:size(dv, 2)
            for kk = 1:size(dv, 3) 
                if dv(ii, jj, kk) < min
                    min = dv(ii, jj, kk);
                    pos1 = ii;
                    pos2 = jj;
                    pos3 = kk;
                end
            end
        end
    end
end

end

