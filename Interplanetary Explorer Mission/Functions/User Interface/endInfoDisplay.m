function [] = endInfoDisplay(data, results)

% Creates the output for the mission
%
% Usage
% [] = endInfoDisplay(data, results)
%
% Input arguments:
% ----------------------------------------------------------------
% data          [-]     data (other function output)        [-]
% results       [-]     results (other function output)     [-]
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

date1 = mjd20002date(results.tspan(1));
date2 = mjd20002date(results.tspan(2));
date3 = mjd20002date(results.tspan(3));

disp("Results: ")
fprintf("Delta V: ")
fprintf("\n - Delta V total: \t\t\t\t\t" + results.dvMin + " km/s")
fprintf("\n - Delta V departure: \t\t\t\t" + results.dv1 + " km/s")
fprintf("\n - Delta V perigee flyby: \t\t\t" + results.dv2 + " km/s")
fprintf("\n - Delta V flyby (given by planet): " + data.deltaVtot + " km/s")
fprintf("\n - Delta V arrival: \t\t\t\t" + results.dv3 + " km/s")

fprintf("\nDelta t: ")
fprintf("\n - Delta t total: \t\t\t\t\t" + (results.tspan(3) - results.tspan(1))/365.25 + " years")
fprintf("\n - Time departure: \t\t\t\t\t" + date1(1) + "/" + date1(2) + "/" + date1(3) + " " + date1(4) + ":" + date1(5) + ":" + date1(6) + " [date]")
fprintf("\n - Time flyby: \t\t\t\t\t\t" + date2(1) + "/" + date2(2) + "/" + date2(3) + " " + date2(4) + ":" + date2(5) + ":" + date2(6) + " [date]")
fprintf("\n - Time arrival: \t\t\t\t\t" + date3(1) + "/" + date3(2) + "/" + date3(3) + " " + date3(4) + ":" + date3(5) + ":" + date3(6) + " [date]")

fprintf("\nFlyby data: ")
fprintf("\n - Perigee radius: \t\t\t\t\t" + data.rp + " km")
fprintf("\n - Planet radius: \t\t\t\t\t" + data.r_planet + " km")
fprintf("\n - Perigee height: \t\t\t\t\t" + data.h + " km")
fprintf("\n - Turn angle: \t\t\t\t\t\t" + rad2deg(data.deltaTot) + " deg")
fprintf("\n - SoI time of flight: \t\t\t\t" + data.ToF/(24*3600) + " days")

fprintf("\n")

end

