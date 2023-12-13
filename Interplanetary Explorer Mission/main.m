%% data

clc, clear
close all

addpath("Functions\")
addpath("Functions\time\")
addpath("Functions_custom\")

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058


%% orbit plot

clc
close all


time_instant = [2028 01 01 0 0 0];
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


%% porkchop plots - departure, fly by

clc, clear
close all

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

departure.planetId = 6;
flyby.plnetId = 5;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);

% dep_time_vect = linspace(date2mjd2000(dep_time), 1.45e4, 100);
% arr_time_vect = linspace(1.47e4, 1.78e4, 100);

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
        [flyby.kep, ~] = uplanet(arr_time_vect(j), flyby.plnetId);

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

dv = dv_1 + dv_2;

contour(dep_time_vect, arr_time_vect, dv', 2:0.2:8)
colorbar, grid on, hold on
xlabel("Departure")
ylabel("Arrival")

[pos1, pos2] = find(dv == min(min(dv)));
plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1)

% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")

% ----------Larmbert arc departure flyby

%Define actual Transfer Mission Parameters

departure.Departure_Date = mjd20002date(ep_time_vect(pos1));
flyby.ArrivalDate = mjd20002date(arr_time_vect(pos2));
ToF_dep_flyby=(ArrTd(arrmin)-DepTd(depmin))*24*60*60;

[departure.kep, ksun] = uplanet(departure.Departure_Date, departure.planetId);
[flyby.kep, ~] = uplanet(flyby.ArrivalDate, flyby.plnetId);


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
    

% %Propagation transfer ARC
% y0 = [ rEreal VIreal ];
% % Set time span
% Tt = 2*pi*sqrt( Areal^3/muSun ); % Orbital period [s]
% tspan0 = linspace( 0,ToFreal, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% % Set options for the ODE solver
% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% % Perform the integration
% [   t, Y ] = ode113( @(t,y) ode_2bp(t,y,muSun), tspan0, y0, options );
% 
% 
% % %Propagation Derparture Range
% y0ED = [ rE vE ];
% % Set time span
% tspanED = linspace( 0,(-dep2d+dep1d)*24*60*60 , 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% % Set options for the ODE solver
% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% % Perform the integration
% [   t, YED ] = ode113( @(t,y) ode_2bp(t,y,ksunE), tspanED, y0ED, options );
% 
% % %Propagation Arrival Range
% y0MA = [ rM vM ];
% % Set time span
% tspanMA = linspace( 0,(-arr2d+arr1d)*24*60*60 , 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% % Set options for the ODE solver
% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% % Perform the integration
% [   t, YMA ] = ode113( @(t,y) ode_2bp(t,y,ksunM), tspanMA, y0MA, options );
% 
% % %Propagation EARTH Motion During Transfer
% y0Emot = [ rEreal vEreal ];
% % Set time span
% tspanEmot = linspace( 0, ToFreal, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% % Set options for the ODE solver
% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% % Perform the integration
% [   t, YEmot ] = ode113( @(t,y) ode_2bp(t,y,ksunEreal), tspanEmot, y0Emot, options );
% 
% % %Propagation MARS orbit
% y0Mmot = [ rMrealfinal vMrealfinal ];
% % Set time span
% TMmot = 2*pi*sqrt( aMreal^3/ksunMreal ); % Orbital period [s]
% tspanMmot = linspace( 0,ToFreal, 5000 ); %% change 2*T to 5*T to get 5 orbital periods
% % Set options for the ODE solver
% options = odeset( 'RelTol', 1e-13, 'AbsTol', 1e-14 );
% % Perform the integration
% [   t, YMmot ] = ode113( @(t,y) ode_2bp(t,y,ksunMreal), tspanMmot, y0Mmot, options );
% 
% 
% % Plot the results
% figure(2)
% plot3( YE(:,1), YE(:,2), YE(:,3), '--','color', 'b' )
% xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
% title('orbits');
% axis equal;
% grid on;
% hold on
% plot3( YM(:,1), YM(:,2), YM(:,3), '--' ,'color', 'r')
% plot3( Y(:,1), Y(:,2), Y(:,3), '-','color', 'g')
% plot3( YED(:,1), YED(:,2), YED(:,3), '-' )
% plot3( YMA(:,1), YMA(:,2), YMA(:,3), '-' )
% plot3( YEmot(:,1), YEmot(:,2), YEmot(:,3), '-','color', 'b' )
% plot3( YMmot(:,1), YMmot(:,2), YMmot(:,3), '-','color', 'r' )
% plot3(rEreal(1),rEreal(2),rEreal(3),'o')
% plot3(rMreal(1),rMreal(2),rMreal(3),'o')
% plot3(0,0,0,'o','color', 'y')
% legend('Earth Orbit','Mars Orbit','Trasnfer Arc','Departure Range','Arrival Range' ...
%     ,'Earth Motion During Transfer','Mars Motion During Transfer')
% hold off
% 






%% porkchop plots - fly by, arrival

clc, clear
close all

dep_time = [2028 01 01 0 0 0];
arr_time = [2058 01 01 0 0 0];

flyby.planetId = 5;
arrival.bodyId = 79;

dep_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);
arr_time_vect = linspace(date2mjd2000(dep_time), date2mjd2000(arr_time), 300);

% dep_time_vect = linspace(1.35e4, 1.48e4, 100);
% arr_time_vect = linspace(1.535e4, 1.565e4, 100);

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

dv = dv_1 + dv_2;

contour(dep_time_vect, arr_time_vect, dv', 12:0.5:25) % 12:35
colorbar, grid on, hold on, axis equal
xlabel("Departure")
ylabel("Arrival")

[pos1, pos2] = find(dv == min(min(dv)));
plot(dep_time_vect(pos1), arr_time_vect(pos2), 'xr', LineWidth=1)

% surface(dep_time_vect, arr_time_vect, dv', EdgeColor="none")



%% Grid Search for departure-flyby-arrival 

clc, clear
close all

%--------Define well this values with Tsyn and ToF!!!!!!
departure.dep_time = [2028 01 01 0 0 0];
departure.arr_time = [2058 01 01 0 0 0];
flyby.arr_time = [2028 01 01 0 0 0];
flyby.dep_time = [2058 01 01 0 0 0];
arrival.dep_time = [2028 01 01 0 0 0];
arrival.arr_time = [2058 01 01 0 0 0];
rp=1000; %[km]
%-----------------------------------------------------

departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;
[flyby.kep, ksun] = uplanet(date2mjd2000(flyby.arr_time), flyby.planetId);
[flyby.r0, flyby.v0] = kep2car([flyby.kep, ksun]);

time1=date2mjd2000(departure.dep_time);
time2=date2mjd2000(departure.arr_time);
time3=date2mjd2000(flyby.dep_time);
time4=date2mjd2000(flyby.arr_time);
time5=date2mjd2000(arrival.dep_time);
time6=date2mjd2000(arrival.arr_time);

%first lambert arc
departure.time_vect = linspace(time1,time2, 20);

flyby.time_vect = linspace(time3,time4, 20);

arrival.time_vect = linspace(time5,time6, 20);


orbitType = 0;
dv_1 = zeros(length(departure.time_vect), length(flyby.time_vect),length(arrival.time_vect));
dv_2 = dv_1;
dv_3 = dv_2;
tof_vect_1 = dv_1;
tof_vect_2 = dv_3;
fixedtol=1;
tol=departure.time_vect(end)-departure.time_vect(1);

while fixedtol<tol
    for i = 1:length(departure.time_vect)
        for j = 1:length(flyby.time_vect)
            for k = 1:length(arrival.time_vect)
    
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
                [flyby.kep, ~] = uplanet(flyby.time_vect(j), flyby.planetId);
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
        
                dv_1(i, j, k) = norm(VI_1 - departure.v0);
                dv_2(i, j, k) = abs(sqrt((2*astroConstants(13)/rp)+norm(VF_1+flyby.v0)^2)-sqrt((2*astroConstants(13)/rp)+norm(VI_2+flyby.v0)^2));
                dv_3(i, j, k) = norm(arrival.v0 - VF_2);
                
            end
        end
    end


dv = dv_1 + dv_2 + dv_3;

[min, pos1, pos2, pos3] = findMin(dv);




if pos1==1
    departure.time_vect=linspace(departure.time_vect(pos1),departure.time_vect(pos1+1),20);
    tol=departure.time_vect(pos1+1)-departure.time_vect(pos1);
else
    departure.time_vect=linspace(departure.time_vect(pos1-1),departure.time_vect(pos1+1),20);
    tol=departure.time_vect(pos1+1)-departure.time_vect(pos1-1);
end

if pos2==1
    flyby.time_vect=linspace(flyby.time_vect(pos2),flyby.time_vect(pos2+1),20);
else
    flyby.time_vect=linspace(flyby.time_vect(pos2-1),flyby.time_vect(pos2+1),20);
end

if pos3==1
    arrival.time_vect=linspace(arrival.time_vect(pos3),arrival.time_vect(pos3+1),20);
else
    arrival.time_vect=linspace(arrival.time_vect(pos3-1),arrival.time_vect(pos3+1),20);

end

disp ("Iteration done!")
end

%% Synodic Periods

clc, clear
close all

time_instant = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(time_instant);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.planetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);


departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [1/s]

flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [1/s]

arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [1/s]

%Synodic Period between departure planet and flyby
T_syn_dep_flyby=(departure.T_orb*flyby.T_orb)/(abs(departure.T_orb-flyby.T_orb)); %[s]

%Synodic Period between flyby planet and arrival
T_syn_flyby_arr=(arrival.T_orb*flyby.T_orb)/(abs(arrival.T_orb-flyby.T_orb)); %[s]

%Synodic Period between departure planet and arrival
T_syn_dep_arr=(arrival.T_orb*departure.T_orb)/(abs(arrival.T_orb-departure.T_orb)); %[s]

%% function find dvmin

function [min, pos1, pos2, pos3] = findMin(dv)
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
                    pos3=k;
                end
            end
        end
    end

    if pos1 == 0 && pos2 == 0 && pos3 == 0
        error("Unable to find minimum")
    end
end
