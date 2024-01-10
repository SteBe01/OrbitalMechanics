function [] = startInfoDisplay(mission)

% Print the data for the chosen mission
%
% Usage
% [] = startInfoDisplay(mission)
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

disp("Mission details: ")
fprintf("Departure Id: " + mission.departure_Id + ", flyby Id: " + mission.flyby_Id + " and Arrival Id: " + mission.arrival_Id)
fprintf("\nWindow details [date]: ")
fprintf("\t - departure: \tfrom " + mission.dep_time_lb(1)+"/"+mission.dep_time_lb(2)+"/"+mission.dep_time_lb(3)+" "+mission.dep_time_lb(4)+":"+mission.dep_time_lb(5)+":"+mission.dep_time_lb(6) + " to " + mission.dep_time_ub(1)+"/"+mission.dep_time_ub(2)+"/"+mission.dep_time_ub(3)+" "+mission.dep_time_ub(4)+":"+mission.dep_time_ub(5)+":"+mission.dep_time_ub(6))
fprintf("\n\t\t\t\t\t\t - flyby: \t\tfrom " + mission.flyby_time_lb(1)+"/"+mission.flyby_time_lb(2)+"/"+mission.flyby_time_lb(3)+" "+mission.flyby_time_lb(4)+":"+mission.flyby_time_lb(5)+":"+mission.flyby_time_lb(6) + " to " + mission.flyby_time_ub(1)+"/"+mission.flyby_time_ub(2)+"/"+mission.flyby_time_ub(3)+" "+mission.flyby_time_ub(4)+":"+mission.flyby_time_ub(5)+":"+mission.flyby_time_ub(6))
fprintf("\n\t\t\t\t\t\t - arrival: \tfrom " + mission.arr_time_lb(1)+"/"+mission.arr_time_lb(2)+"/"+mission.arr_time_lb(3)+" "+mission.arr_time_lb(4)+":"+mission.arr_time_lb(5)+":"+mission.arr_time_lb(6) + " to " + mission.arr_time_ub(1)+"/"+mission.arr_time_ub(2)+"/"+mission.arr_time_ub(3)+" "+mission.arr_time_ub(4)+":"+mission.arr_time_ub(5)+":"+mission.arr_time_ub(6))
fprintf("\n")

end

