function [] = windowDisplay(mission)

% Print the time window for the chosen mission
%
% Usage
% [] = windowDisplay(mission)
%
% Input arguments:
% ----------------------------------------------------------------
% mission           [-]     mission struct      [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% N/A
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

mission.dep_time_lb = date2mjd2000(mission.dep_time_lb);
mission.dep_time_ub = date2mjd2000(mission.dep_time_ub);
mission.flyby_time_lb = date2mjd2000(mission.flyby_time_lb);
mission.flyby_time_ub = date2mjd2000(mission.flyby_time_ub);
mission.arr_time_lb = date2mjd2000(mission.arr_time_lb);
mission.arr_time_ub = date2mjd2000(mission.arr_time_ub);

porkchop_start(mission.dep_time, mission.arr_time, mission.departure_Id, mission.flyby_Id, mission.arrival_Id);
sgtitle("Chosen search window")

subplot(1, 2, 1)
plot([mission.dep_time_lb mission.dep_time_ub], [mission.flyby_time_lb mission.flyby_time_lb], 'r');
hold on
plot([mission.dep_time_lb mission.dep_time_ub], [mission.flyby_time_ub mission.flyby_time_ub], 'r');
plot([mission.dep_time_lb mission.dep_time_lb], [mission.flyby_time_lb mission.flyby_time_ub], 'r');
plot([mission.dep_time_ub mission.dep_time_ub], [mission.flyby_time_lb mission.flyby_time_ub], 'r');
subplot(1, 2, 2)
plot([mission.flyby_time_lb mission.flyby_time_ub], [mission.arr_time_lb mission.arr_time_lb], 'r');
hold on
plot([mission.flyby_time_lb mission.flyby_time_ub], [mission.arr_time_ub mission.arr_time_ub], 'r');
plot([mission.flyby_time_lb mission.flyby_time_lb], [mission.arr_time_lb mission.arr_time_ub], 'r');
plot([mission.flyby_time_ub mission.flyby_time_ub], [mission.arr_time_lb mission.arr_time_ub], 'r');

end

