%% Grid Search with smaller windows

% This script was used to create the "Grid Search", later implemented
% as a function and moved to the main script

clc, clear
close all

windowType = 1;

restoredefaultpath
addpath(genpath("..\"))

xcust(1) = 1.289338494887906e+04;
xcust(2) = 1.670686199971615e+04;
xcust(3) = 1.754357184719549e+04;

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

[departureTime, flybyTime, arrivalTime] = loadWindows(windowType);

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

window_size = 30;
fixedtol = 1e3;                         % in seconds
fmincon_choice = 3;                     % 0 for no fmincon
orbitType = 0;
animation = 1;

if animation
    porkchop_start(dep_time, arr_time, departure.planetId, flyby.planetId, arrival.bodyId);
end

for totWindows = 1:length(departureTime)
    departure.time_vect = linspace(departureTime{totWindows}.lb, departureTime{totWindows}.ub, window_size);
    flyby.time_vect = linspace(flybyTime{totWindows}.lb, flybyTime{totWindows}.ub, window_size);
    arrival.time_vect = linspace(arrivalTime{totWindows}.lb, arrivalTime{totWindows}.ub, window_size);
    
    dv_1 = NaN* ones(length(departure.time_vect), length(flyby.time_vect), length(arrival.time_vect));
    dv_2 = dv_1;
    dv_3 = dv_1;
    rp = dv_1;
    num_iter = 0;
    
    tol = (departure.time_vect(end)-departure.time_vect(1)) * 24 * 60 * 60;
    
    while fixedtol < tol
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
            
                    [dv_1(i, j, k), dv_2(i, j, k), dv_3(i, j, k), rp_temp, exitValue, ~] = completeInterplanetary(departure.time_vect(i), flyby.time_vect(j), arrival.time_vect(k), departure.planetId, flyby.planetId, arrival.bodyId);
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
    
        if tol > fixedtol
            departure.time_vect = linspace(departure.time_vect(pos1_bnd_1), departure.time_vect(pos1_bnd_2), window_size);
            flyby.time_vect = linspace(flyby.time_vect(pos2_bnd_1), flyby.time_vect(pos2_bnd_2), window_size);
            arrival.time_vect = linspace(arrival.time_vect(pos3_bnd_1), arrival.time_vect(pos3_bnd_2), window_size);
        end
    
        if fmincon_choice == num_iter
            break
        else
            disp ("iteration " + num_iter + " done in " + toc + " s, dv = " + dVmin + " km/s")
        end
    end
    
    if fmincon_choice ~= 0
        opts = optimset('TolX', eps(1), 'TolFun', eps(1), 'Display', 'off');
        lb = date2mjd2000(dep_time) * ones(3,1);
        ub = date2mjd2000(arr_time) * ones(3,1);
        [tspan, dv_fmin] = fmincon(@(tspan) completeInterplanetaryGS(tspan(1), tspan(2), tspan(3), departure.planetId, flyby.planetId, arrival.bodyId), [departure.time_vect(pos1) flyby.time_vect(pos2) arrival.time_vect(pos3)]', [], [], [], [], lb, ub, [], opts);
        
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
    a1 = plot([departureTime{totWindows}.lb departureTime{totWindows}.ub], [flybyTime{totWindows}.lb flybyTime{totWindows}.lb], 'r', LineWidth=1);
    hold on
    a2 = plot([departureTime{totWindows}.lb departureTime{totWindows}.ub], [flybyTime{totWindows}.ub flybyTime{totWindows}.ub], 'r', LineWidth=1);
    a3 = plot([departureTime{totWindows}.lb departureTime{totWindows}.lb], [flybyTime{totWindows}.lb flybyTime{totWindows}.ub], 'r', LineWidth=1);
    a4 = plot([departureTime{totWindows}.ub departureTime{totWindows}.ub], [flybyTime{totWindows}.lb flybyTime{totWindows}.ub], 'r', LineWidth=1);
    a5 = plot(tspan(1), tspan(2), 'xr', LineWidth=1);
    subplot(1, 2, 2)
    b1 = plot([flybyTime{totWindows}.lb flybyTime{totWindows}.ub], [arrivalTime{totWindows}.lb arrivalTime{totWindows}.lb], 'r', LineWidth=1);
    hold on
    b2 = plot([flybyTime{totWindows}.lb flybyTime{totWindows}.ub], [arrivalTime{totWindows}.ub arrivalTime{totWindows}.ub], 'r', LineWidth=1);
    b3 = plot([flybyTime{totWindows}.lb flybyTime{totWindows}.lb], [arrivalTime{totWindows}.lb arrivalTime{totWindows}.ub], 'r', LineWidth=1);
    b4 = plot([flybyTime{totWindows}.ub flybyTime{totWindows}.ub], [arrivalTime{totWindows}.lb arrivalTime{totWindows}.ub], 'r', LineWidth=1);
    b5 = plot(tspan(2), tspan(3), 'xr', LineWidth=1);
    sgtitle("Delta velocity = " + dv_fmin + " km/s")
    drawnow

    solutions.dvMin{totWindows} = dv_fmin;
    solutions.tspan{totWindows} = tspan;
end


%% Plot

n = 14;
if windowType == 1
    n = n - 12;
end

departure_time = solutions.tspan{n}(1);
flyby_time = solutions.tspan{n}(2);
arrival_time = solutions.tspan{n}(3);

missionPlot(departure_time, flyby_time, arrival_time, departure.planetId, flyby.planetId, arrival.bodyId)
data = flybyPlot(departure_time, flyby_time, arrival_time, departure.planetId, flyby.planetId, arrival.bodyId, 1e7);


%% Functions

function [min, pos1, pos2, pos3] = findMin3(dv)
    min = max(max(max(dv)));
    pos1 = 0;
    pos2 = 0;
    pos3 = 0;
    for i = 1:size(dv, 1) 
        for j = 1:size(dv, 2)
            for k = 1:size(dv, 3) 
                if dv(i, j, k) < min
                    min = dv(i, j, k);
                    pos1 = i;
                    pos2 = j;
                    pos3 = k;
                end
            end
        end
    end
end

