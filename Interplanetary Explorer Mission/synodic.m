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

