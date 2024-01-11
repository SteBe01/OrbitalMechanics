function [a_sat] = aFinder(k, m, om_E, mu)

% a for repeating orbits
%
% Usage
% [a_sat] = aFinder(k, m, om_E, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% k             [1x1]   satellite revolutions       [-]
% m             [1x1]   Earth revolutions           [-]
% om_E          [1x1]   Earth angular velocity      [deg/h]
% mu            [1x1]   const                       [km^3/s^2]
% 
% Output arguments:
% -----------------------------------------------------------------
% a             [1x1]   semi-major                  [km]
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

n = (deg2rad(om_E)/3600) * (k / m);
a_sat = (mu/n^2)^(1/3);

end
