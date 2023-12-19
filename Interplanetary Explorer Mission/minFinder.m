%% data

clc, clear
close all

restoredefaultpath
addpath(genpath("Functions\"))
addpath("Functions_custom\")

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058


%% combined porkchops

clc, clear
close all

warning("Not working!!!!!!!!!!!!!!!!, accrocchio")

xcust(1) = 1.554407711683449e+04;
xcust(2) = 1.878584012172524e+04;
xcust(3) = 1.968280293314873e+04;

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);

orbitType = 0;
dv_1 = zeros(length(dep_time_vect), length(arr_time_vect));
for i = 1:length(dep_time_vect)
    for j = i:length(arr_time_vect)
        tof = (arr_time_vect(j) - dep_time_vect(i)) * 24 * 60 * 60; % seconds

        if tof <= 1e4
            dv_1(i, j) = NaN;
            continue
        end

        [departure.kep, ksun] = uplanet(dep_time_vect(i), departure.planetId);
        [flyby.kep, ~] = uplanet(arr_time_vect(j), flyby.planetId);

        [departure.r0, departure.v0] = kep2car([departure.kep ksun]);
        [flyby.r0, ~] = kep2car([flyby.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(departure.r0, flyby.r0, tof, ksun, orbitType, 0);
        if A < 0
            dv_1(i, j) = NaN;
            continue
        end

        dv_1(i, j) = norm(VI - departure.v0);
    end
end

dv_2 = zeros(length(dep_time_vect), length(arr_time_vect));
for i = 1:length(dep_time_vect)
    for j = i:length(arr_time_vect)
        tof = (arr_time_vect(j) - dep_time_vect(i)) * 24 * 60 * 60; % seconds

        if tof <= 1e4
            dv_2(i, j) = NaN;
            continue
        end

        [flyby.kep, ksun] = uplanet(dep_time_vect(i), flyby.planetId);
        [arrival.kep, ~] = ephNEO(arr_time_vect(j), arrival.bodyId);

        [flyby.r0, ~] = kep2car([flyby.kep ksun]);
        [arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(flyby.r0, arrival.r0, tof, ksun, orbitType, 0);
        if A < 0
            dv_2(i, j) = NaN;
            continue
        end

        dv_2(i, j) = norm(arrival.v0 - VF);
    end
end


% merge matrix - rip
test1 = NaN * ones(length(dv_1), length(dv_2));
k = 1;
for i = 1:length(dv_1)
    for j = i:length(dv_2)
        dv_matrix_dv1 = dv_1(i, j);

        if isnan(dv_matrix_dv1)
            continue
        end

        [min, pos1, pos2, flag] = findMin2(dv_2(:, j));
        if ~flag
            continue
        end

        [~, dv21, ~, rp, exitValue] = completeInterplanetary(dep_time_vect(i), dep_time_vect(pos2), dep_time_vect(j), departure.planetId, flyby.planetId, arrival.bodyId);
        if dv21 > 0.6
            continue
        end

        num.dv_matrix_dv1{k} = dv_matrix_dv1 + min;
        num.index_dep{k} = dep_time_vect(i);
        num.index_flyby{k} = dep_time_vect(j);
        num.index_arrival{k} = dep_time_vect(pos1);
        num.index1{k} = i;
        num.index2{k} = j;
        num.index3{k} = pos2;
        k = k + 1;

        test1(300-i, pos1) = dv_matrix_dv1 + min;
    end
end

contour(dep_time_vect, arr_time_vect, mirror(test1), 1:0.1:15,HandleVisibility="off") % 12:35
colorbar, grid on, hold on, axis equal
plot(xcust(1), xcust(3), 'or', LineWidth=1)
xlabel("Departure")
ylabel("Arrival")



%% combined porkchops - aaaaaaaaaaa

clc, clear
close all

xcust(1) = 1.554407711683449e+04;
xcust(2) = 1.878584012172524e+04;
xcust(3) = 1.968280293314873e+04;

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 200);

orbitType = 0;
dv_1 = NaN .* ones(length(time_vect), length(time_vect));
test1 = dv_1;
for i = 1:length(time_vect)
    for j = i:length(time_vect)
        tof = (time_vect(j) - time_vect(i)) * 24 * 60 * 60; % seconds

        if tof < 1e3
            continue
        end

        [departure.kep, ksun] = uplanet(time_vect(i), departure.planetId);
        [flyby.kep, ~] = uplanet(time_vect(j), flyby.planetId);

        [departure.r0, departure.v0] = kep2car([departure.kep ksun]);
        [flyby.r0, ~] = kep2car([flyby.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(departure.r0, flyby.r0, tof, ksun, orbitType, 0);
        if A < 0
            continue
        end

        dv_1(i, j) = norm(VI - departure.v0);

        dv_2 = NaN .* ones(length(time_vect), length(time_vect));
        for l = j:length(time_vect)
            tof1 = (time_vect(l) - time_vect(j)) * 24 * 60 * 60; % seconds
    
            [flyby.kep, ksun] = uplanet(time_vect(j), flyby.planetId);
            [arrival.kep, ~] = ephNEO(time_vect(l), arrival.bodyId);
    
            [flyby.r0, ~] = kep2car([flyby.kep ksun]);
            [arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);
    
            [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(flyby.r0, arrival.r0, tof1, ksun, orbitType, 0);
            if A < 0
                continue
            end
    
            dv_2(j, l) = norm(arrival.v0 - VF);
        end
        
        [dv1, dv2, dv3, rp, exitValue] = completeInterplanetary(time_vect(i), time_vect(j), time_vect(l), departure.planetId, flyby.planetId, arrival.bodyId);

        if dv2 > 2
            continue
        end
        [min, pos1, pos2, flag] = findMin2(dv_2);
        if flag && ~isnan(dv_1(i,j))
            % disp(pos2 + " " + i)
            test1(i, pos1) = min + dv_1(i, j);

        end
    end
end
%%
contour(time_vect, time_vect, test1', 10:0.1:18, HandleVisibility="off") % 12:35
colorbar, grid on, hold on, axis equal
plot(xcust(1), xcust(3), 'or', LineWidth=1)
xlabel("Departure")
ylabel("Arrival")



%% animation

clc, clear
close all

%--------Define well this values with Tsyn and ToF!!!!!!
mission.dep_time = [2028 01 01 0 0 0];
mission.arr_time = [2058 01 01 0 0 0];
%-----------------------------------------------------

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;
fixedtol = 1e3;
window_size = 50;

time1 = date2mjd2000(mission.dep_time);
time2 = date2mjd2000(mission.arr_time);

%first lambert arc
departure.time_vect = linspace(time1, time2, window_size);
flyby.time_vect = linspace(time1, time2, window_size);
arrival.time_vect = linspace(time1, time2, window_size);

orbitType = 0;
dv_1 = zeros(length(departure.time_vect), length(flyby.time_vect), length(arrival.time_vect));
dv_2 = dv_1;
dv_3 = dv_1;
dv_flyby_tot = dv_1;
tof_vect_1 = dv_1;
tof_vect_2 = dv_1;
num_iter = 0;

tol = (departure.time_vect(end)-departure.time_vect(1)) * 24 * 60 * 60;


num_iter = num_iter + 1;
iteration.number{num_iter} = num_iter;

for i = 1:length(departure.time_vect)
    for j = i:length(flyby.time_vect)
        for k = j:length(arrival.time_vect)

            tof_1 = (flyby.time_vect(j) - departure.time_vect(i)) * 24 * 60 * 60; % seconds
            tof_2 = (arrival.time_vect(k) - flyby.time_vect(j)) * 24 * 60 * 60; % seconds
    
            if (tof_1 <= 1e5) && (tof_2 <= 1e5)
                dv_1(i, j, k) = NaN;
                dv_2(i, j, k) = NaN;
                dv_3(i, j, k) = NaN;
                continue
            end
            tof_vect_1(i, j, k) = tof_1;
            tof_vect_2(i, j, k) = tof_2;
    
            [departure.kep, ksun] = uplanet(departure.time_vect(i), departure.planetId);
            [flyby.kep, flyby.mu] = uplanet(flyby.time_vect(j), flyby.planetId);
            [arrival.kep, ~] = ephNEO(arrival.time_vect(k), arrival.bodyId);
    
            [departure.r0, departure.v0] = kep2car([departure.kep ksun]);
            [flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);
            [arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);
    
            [A_1, P_1, E_1, ERROR_1, VI_1, VF_1, TPAR_1, THETA_1] = lambertMR(departure.r0, flyby.r0, tof_1, ksun, orbitType, 0);
            if A_1 < 0
                dv_1(i, j, k) = NaN;
                dv_2(i, j, k) = NaN;
                dv_3(i, j, k) = NaN;
                continue
            end

            [A_2, P_2, E_2, ERROR_2, VI_2, VF_2, TPAR_2, THETA_2] = lambertMR(flyby.r0, arrival.r0, tof_2, ksun, orbitType, 0);
            if A_2 < 0
                dv_1(i, j, k) = NaN;
                dv_2(i, j, k) = NaN;
                dv_3(i, j, k) = NaN;
                continue
            end
    
            v_inf_minus = VI_2 + flyby.v0;
            v_inf_plus = VF_1 + flyby.v0;
            dv_1(i, j, k) = norm(VI_1 - departure.v0);
            [rp, flag] = rpsolver(v_inf_minus, v_inf_plus, flyby.planetId);
            if flag == 0
                dv_1(i, j, k) = NaN;
                dv_2(i, j, k) = NaN;
                dv_3(i, j, k) = NaN;
                continue
            end
            dv_2(i, j, k) = abs(sqrt((2*astroConstants(flyby.planetId + 10)/rp)+norm(v_inf_plus)^2)-sqrt((2*astroConstants(flyby.planetId + 10)/rp)+norm(v_inf_minus)^2));
            dv_3(i, j, k) = norm(arrival.v0 - VF_2);
            dv_flyby_tot(i, j, k) = norm(VF_1 - VF_2);
        end
    end
end

dv = dv_1 + dv_2 + dv_3;

figure

for i = 1:50
    array2plot = reshape(dv(:,i,:), window_size, window_size);
    a = surf(departure.time_vect, arrival.time_vect, array2plot');
    colorbar, grid on, hold on
    xlabel("Departure")
    ylabel("Arrival")
    drawnow
    view(340, 50)
    pause(1)
    delete(a)
    disp(i)
end



%% function find dvmin

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

    if pos1 == 0 && pos2 == 0 && pos3 == 0
        error("Unable to find minimum")
    end
end


function [min, pos1, pos2, flag] = findMin2(dv)
    min = max(max(max(dv)));
    pos1 = 0;
    pos2 = 0;
    for i = 1:size(dv, 1) 
        for j = 1:size(dv, 2)
            if dv(i, j) < min
                min = dv(i, j);
                pos1 = i;
                pos2 = j;
            end
        end
    end

    flag = 1;
    if pos1 == 0 && pos2 == 0
        flag = 0;
    end
end


function [array] = mirror(array)
    temp = array;
    for i = 1:length(temp)
        array(i, :) = temp(length(temp) - i + 1, :);
    end
    for i = 1:length(temp)
        array(:, i) = temp(:, length(temp) - i + 1);
    end
end

