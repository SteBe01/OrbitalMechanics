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


%% orbits plot

clc
close all

time_instant = [2028 01 01 0 0 0];          % TBD
departure.planetId = 6;
flyby.plnetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(time_instant);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.plnetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);


% Set options for the ODE solver
options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );

% orbit propagation - departure
[departure.r0, departure.v0] = kep2car([departure.kep, ksun]);
departure.y0 = [departure.r0 departure.v0];

departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [1/s]
tspan= linspace( 0, departure.T_orb, 200 );
[ ~, departure.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, departure.y0, options );

% orbit propagation - flyby
[flyby.r0, flyby.v0] = kep2car([flyby.kep, ksun]);
flyby.y0 = [flyby.r0 flyby.v0];

flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [1/s]
tspan= linspace( 0, flyby.T_orb, 200 );
[ ~, flyby.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, flyby.y0, options );

% orbit propagation - arrival
[arrival.r0, arrival.v0] = kep2car([arrival.kep, ksun]);
arrival.y0 = [arrival.r0 arrival.v0];

arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [1/s]
tspan= linspace( 0, arrival.T_orb, 200 );
[ ~, arrival.Y ] = ode113( @(t,y) ode_2bp(t,y,ksun), tspan, arrival.y0, options );


% plots
figure
axis equal, grid on, hold on
plot3(departure.Y(:, 1), departure.Y(:, 2), departure.Y(:, 3), Color="blu", HandleVisibility="off")
plot3(flyby.Y(:, 1), flyby.Y(:, 2), flyby.Y(:, 3), Color="green", HandleVisibility="off")
plot3(arrival.Y(:, 1), arrival.Y(:, 2), arrival.Y(:, 3), Color="red", HandleVisibility="off")

scatter3(departure.r0(1), departure.r0(2), departure.r0(3), "blu", "filled")
scatter3(flyby.r0(1), flyby.r0(2), flyby.r0(3), "green", "filled")
scatter3(arrival.r0(1), arrival.r0(2), arrival.r0(3), "red", "filled")


legend("Departure (Saturn)", "Fly-By (Jupiter)", "Arrival (Asteroid N.79)", Location="best")

xlim([min(min(min(departure.Y(:,1)), min(flyby.Y(:,1))), min(arrival.Y(:,1))) max(max(max(departure.Y(:,1)), max(flyby.Y(:,1))), max(arrival.Y(:,1)))])
ylim([min(min(min(departure.Y(:,2)), min(flyby.Y(:,2))), min(arrival.Y(:,2))) max(max(max(departure.Y(:,2)), max(flyby.Y(:,2))), max(arrival.Y(:,2)))])
zlim([min(min(min(departure.Y(:,3)), min(flyby.Y(:,3))), min(arrival.Y(:,3))) max(max(max(departure.Y(:,3)), max(flyby.Y(:,3))), max(arrival.Y(:,3)))])

xlabel("x [km]")
ylabel("y [km]")
zlabel("z [km]")

view(30,30)



%% porkchop plots - departure - fly by and flyby - arrival

clc, clear
close all

xcust(1) = 1.289338494887906e+04;
xcust(2) = 1.670686199971615e+04;
xcust(3) = 1.754357184719549e+04;

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 500);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 500);

orbitType = 0;
dv_1 = zeros(length(dep_time_vect), length(arr_time_vect));
dv_2 = dv_1;
tof_vect = dv_1;
for i = 1:length(dep_time_vect)
    for j = 1:length(arr_time_vect)
        tof = (arr_time_vect(j) - dep_time_vect(i)) * 24 * 60 * 60; % seconds

        if tof <= 1e5
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end
        tof_vect(i, j) = tof;

        [departure.kep, ksun] = uplanet(dep_time_vect(i), departure.planetId);
        [flyby.kep, ~] = uplanet(arr_time_vect(j), flyby.planetId);

        [departure.r0, departure.v0] = kep2car([departure.kep ksun]);
        [flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(departure.r0, flyby.r0, tof, ksun, orbitType, 0);
        if A < 0
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end

        dv_1(i, j) = norm(VI - departure.v0);
        dv_2(i, j) = norm(flyby.v0 - VF);
    end
end

% dv = dv_1 + dv_2;
dv = dv_1; % without flyby dv

contour(dep_time_vect, arr_time_vect, dv', 2:0.2:8, HandleVisibility="off")
colorbar, grid on, hold on
xlabel("Departure")
ylabel("FlyBy")

[pos1, pos2] = find(dv == min(min(dv)));
plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1)
plot(xcust(1), xcust(2), 'or', LineWidth=1)

legend("Global min", "Mission min", Location="best")

% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")



% porkchop plots - fly by, arrival

orbitType = 0;
dv_1 = zeros(length(dep_time_vect), length(arr_time_vect));
dv_2 = dv_1;
tof_vect = dv_1;
for i = 1:length(dep_time_vect)
    for j = 1:length(arr_time_vect)
        tof = (arr_time_vect(j) - dep_time_vect(i)) * 24 * 60 * 60; % seconds

        if tof <= 1e5
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end
        tof_vect(i, j) = tof;

        [flyby.kep, ksun] = uplanet(dep_time_vect(i), flyby.planetId);
        [arrival.kep, ~] = ephNEO(arr_time_vect(j), arrival.bodyId);

        [flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);
        [arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);

        [A, P, E, ERROR, VI, VF, TPAR, THETA] = lambertMR(flyby.r0, arrival.r0, tof, ksun, orbitType, 0);
        if A < 0
            dv_1(i, j) = NaN;
            dv_2(i, j) = NaN;
            continue
        end

        dv_1(i, j) = norm(VI - flyby.v0);
        dv_2(i, j) = norm(arrival.v0 - VF);
    end
end

% dv = dv_1 + dv_2;
dv = dv_2; % without flyby dv

figure
contour(dep_time_vect, arr_time_vect, dv', 4:0.5:25, HandleVisibility="off") % 12:35
colorbar, grid on, hold on, axis equal
xlabel("FlyBy")
ylabel("Arrival")

[pos1, pos2] = find(dv == min(min(dv)));
plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1)
plot(xcust(2), xcust(3), 'or', LineWidth=1)

legend("Global min", "Mission min", Location="best")

% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")

warning("To do: add same colorbar for the plots")

