%% data

clc, clear
close all

restoredefaultpath
addpath(genpath("Functions\"))
addpath(genpath("Functions_custom\"))

simChoice = 1;

% data - goup 2346
% Departure: Saturn
% Flyby: Jupiter
% Arrival: Asteriod N.79
% Earliest Departure: 00:00:00 01/01/2028
% Latest Arrival: 00:00:00 01/01/2058

mission.departure_Id = 6;
mission.flyby_Id = 5;
mission.arrival_Id = 79;

mission.dep_time_lb = [2030 12 09 0 0 0];
mission.dep_time_ub = [2039 02 25 0 0 0];
mission.flyby_time_lb = [2042 12 26 0 0 0];
mission.flyby_time_ub = [2048 09 25 0 0 0];
mission.arr_time_lb = [2045 06 13 0 0 0];
mission.arr_time_ub = [2058 01 01 0 0 0];

% mission.dep_time_lb = [2028 01 01 0 0 0];
% mission.dep_time_ub = [2058 01 01 0 0 0];
% mission.flyby_time_lb = [2028 01 01 0 0 0];
% mission.flyby_time_ub = [2058 01 01 0 0 0];
% mission.arr_time_lb = [2028 01 01 0 0 0];
% mission.arr_time_ub = [2058 01 01 0 0 0];

if simChoice == 1
    % addpath(genpath("Grid Search\"))

    % grid search options
    mission.options.windowType = 1;
    mission.options.window_size = 30;
    mission.options.fixedtol = 1e3;                         % in seconds
    mission.options.fmincon_choice = 3;                     % 0 for no fmincon
    mission.options.animation = 1;

    [solutions] = gridSearch_function(mission);
elseif simChoice == 2
    % addpath("Genetic Algorithm\Functions\")

    % ga options
    mission.options.n_iter = 1;

    [solutions] = ga_function(mission);
elseif simChoice == 3
    % addpath("Global Search\Functions\")

    % multi start options
    mission.options.n_elements = 1e4;
    mission.options.parallel = 1;

    [solutions] = multiStart_function(mission);
else
    error("Invalid simulation choice")
end











%% Find Departure and Arrival Dates Ranges

% ---Synodic Periods-----------

date_departure = [2028 01 01 0 0 0];
departure.planetId = 6;
flyby.planetId = 5;
arrival.bodyId = 79;

time_instant_mjd200 = date2mjd2000(date_departure);

[departure.kep, ksun] = uplanet(time_instant_mjd200, departure.planetId);
[flyby.kep, ~] = uplanet(time_instant_mjd200, flyby.planetId);
[arrival.kep, ~, ~] = ephNEO(time_instant_mjd200, arrival.bodyId);


departure.T_orb = 2*pi*sqrt( departure.kep(1)^3/ksun ); % Orbital period [1/s]

flyby.T_orb = 2*pi*sqrt( flyby.kep(1)^3/ksun ); % Orbital period [1/s]

arrival.T_orb = 2*pi*sqrt( arrival.kep(1)^3/ksun ); % Orbital period [1/s]

%Synodic Period between departure planet and flyby
T_syn_dep_flyby=(departure.T_orb*flyby.T_orb)/(abs(departure.T_orb-flyby.T_orb)); %[s]
T_syn_dep_flyby_years=T_syn_dep_flyby/(60*60*24*astroConstants(32)); %[years]

%Synodic Period between flyby planet and arrival
T_syn_flyby_arr=(arrival.T_orb*flyby.T_orb)/(abs(arrival.T_orb-flyby.T_orb)); %[s]
T_syn_flyby_arr_years=T_syn_flyby_arr/(60*60*24*astroConstants(32)); %[years]

%Synodic Period between departure planet and arrival
T_syn_dep_arr=(arrival.T_orb*departure.T_orb)/(abs(arrival.T_orb-departure.T_orb)); %[s]
T_syn_dep_arr_years=T_syn_dep_arr/(60*60*24*astroConstants(32)); %[years]

%Global Period 
T_Global=1/(abs((1/T_syn_dep_flyby)-(1/T_syn_flyby_arr)));


%Set possible dates for departure and arrival
% jd_arrival_flyby_planet=mjd20002jd(time_instant_mjd200)+T_syn_dep_flyby/(60*60*24);
% date_flyby = jd2date(jd_arrival_flyby_planet);
% 
% jd_arrival_arrival_planet=jd_arrival_flyby_planet+T_syn_flyby_arr/(60*60*24);
% date_arrival= jd2date(jd_arrival_arrival_planet);
jd_arrival=mjd20002jd(time_instant_mjd200)+T_Global/(60*60*24);
date_arrival= jd2date(jd_arrival);

