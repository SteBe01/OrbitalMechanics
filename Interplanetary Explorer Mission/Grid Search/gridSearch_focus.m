%% Grid Search with lookup

clc, clear
close all

restoredefaultpath
addpath("Functions\")
addpath("..\Functions_custom\")
addpath(genpath("..\\Functions"))

xcust(1) = 1.289338494887906e+04;
xcust(2) = 1.670686199971615e+04;
xcust(3) = 1.754357184719549e+04;

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

% whole window
% departure.lb = 1.15e4;
% departure.ub = 1.42e4;
% flyby.lb = 1.49e4;
% flyby.ub = 1.74e4;
% arrival.lb = 1.52e4;
% arrival.ub = date2mjd2000(arr_time);

[departureTime, flybyTime, arrivalTime] = loadWindows;

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
            
                    [dv_1(i, j, k), dv_2(i, j, k), dv_3(i, j, k), rp_temp, exitValue] = completeInterplanetary(departure.time_vect(i), flyby.time_vect(j), arrival.time_vect(k), departure.planetId, flyby.planetId, arrival.bodyId);
                    if exitValue == 0
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
        A_fmin = [-1 1 0; 0 -1 1]; b_fmin = [0 0];
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

    solutions.dvMin{totWindows} = dv_fmin;
    solutions.tspan{totWindows} = tspan;
end


%% Plot

n = 2;

departure.Date = solutions.tspan{n}(1);
flyby.Date = solutions.tspan{n}(2);
arrival.Date = solutions.tspan{n}(3);

ToF_dep_flyby=(flyby.Date - departure.Date)*24*60*60;
ToF_flyby_arr=(arrival.Date - flyby.Date)*24*60*60;

[departure.kep_actual, ksun_actual] = uplanet(departure.Date, departure.planetId);
[flyby.kep_actual, ~] = uplanet(flyby.Date, flyby.planetId);
[arrival.kep_actual, ~] = ephNEO(arrival.Date, arrival.bodyId);

% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% orbit propagation - departure
[departure.r0_actual, departure.v0_actual] = kep2car([departure.kep_actual, ksun_actual]);
departure.y0_actual = [departure.r0_actual departure.v0_actual];

departure.T_orb_actual = 2*pi*sqrt( departure.kep_actual(1)^3/ksun_actual ); % Orbital period [1/s]
departure.tspan= linspace( 0, departure.T_orb_actual, 200 );
[ ~, departure.Y_actual ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), departure.tspan, departure.y0_actual, options );

% orbit propagation - flyby
[flyby.r0_actual, flyby.v0_actual] = kep2car([flyby.kep_actual, ksun_actual]);
flyby.y0_actual = [flyby.r0_actual flyby.v0_actual];

flyby.T_orb_actual = 2*pi*sqrt( flyby.kep_actual(1)^3/ksun_actual ); % Orbital period [1/s]
flyby.tspan= linspace( 0, flyby.T_orb_actual, 200 );
[ ~, flyby.Y_actual ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), flyby.tspan, flyby.y0_actual, options );

% orbit propagation - arrival
[arrival.r0_actual, arrival.v0_actual] = kep2car([arrival.kep_actual, ksun_actual]);
arrival.y0_actual = [arrival.r0_actual arrival.v0_actual];

arrival.T_orb_actual = 2*pi*sqrt( arrival.kep_actual(1)^3/ksun_actual ); % Orbital period [1/s]
arrival.tspan= linspace( 0, arrival.T_orb_actual, 200 );
[ ~, arrival.Y_actual ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), arrival.tspan, arrival.y0_actual, options );

%Propagation first transfer ARC
[A_1_actual, P_1_actual, E_1_actual, ERROR_1_actual, VI_1_actual, VF_1_actual, TPAR_1_actual, THETA_1_actual] = lambertMR(departure.r0_actual, flyby.r0_actual, ToF_dep_flyby, ksun_actual, orbitType, 0);
y0_1 = [ departure.r0_actual VI_1_actual ];
% Set time span
T_1 = 2*pi*sqrt( A_1_actual^3/ksun_actual ); % Orbital period [s]
tspan_1 = linspace( 0,ToF_dep_flyby, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% Perform the integration
[   t, Y_1 ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), tspan_1, y0_1, options );

%Propagation second transfer ARC
[A_2_actual, P_2_actual, E_2_actual, ERROR_2_actual, VI_2_actual, VF_2_actual, TPAR_2_actual, THETA_2_actual] = lambertMR(flyby.r0_actual, arrival.r0_actual, ToF_flyby_arr, ksun_actual, orbitType, 0);
y0_2 = [ flyby.r0_actual VI_2_actual ];
% Set time span
T_2 = 2*pi*sqrt( A_2_actual^3/ksun_actual ); % Orbital period [s]
tspan_2 = linspace( 0,ToF_flyby_arr, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% Perform the integration
[   t, Y_2 ] = ode113( @(t,y) ode_2bp(t,y,ksun_actual), tspan_2, y0_2, options );

% Plot the results
figure()
plot3( departure.Y_actual(:,1), departure.Y_actual(:,2), departure.Y_actual(:,3), '-','color', 'b' )
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('orbits');
axis equal;
grid on;
hold on
plot3( flyby.Y_actual(:,1),flyby.Y_actual(:,2), flyby.Y_actual(:,3), '-' ,'color', 'r')
plot3( arrival.Y_actual(:,1), arrival.Y_actual(:,2), arrival.Y_actual(:,3), '-','color', 'g')
plot3( Y_1(:,1), Y_1(:,2), Y_1(:,3), '--','color', 'm' )
plot3( Y_2(:,1), Y_2(:,2), Y_2(:,3), '--','color', 'm' )
plot3(departure.r0_actual(1),departure.r0_actual(2),departure.r0_actual(3),'o','Color','b','MarkerFaceColor','b')
plot3(flyby.r0_actual(1),flyby.r0_actual(2),flyby.r0_actual(3),'o','Color','r','MarkerFaceColor','r')
plot3(arrival.r0_actual(1),arrival.r0_actual(2),arrival.r0_actual(3),'o','Color','g','MarkerFaceColor','g')
plot3(0,0,0,'o','Color','y','MarkerFaceColor','y')
legend('Saturn Orbit','Jupiter Orbit','Asteroid N.79 Orbit','Transfer Arc 1','Transfer Arc 2', ...
    'Saturn','Jupiter','Asteroid N.79','Sun')
hold off


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
