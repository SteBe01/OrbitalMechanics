function [deltaT] = tof_soi(a, e, rsoi, mu)

% Computes delta time for 1 hyperbolic leg (between pericenter and rsoi)
%
% Usage
% [deltaT] = tof_soi(a, e, rsoi, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% a             [1x1]       semi-major axis                 [km]
% e             [1x1]       eccentricity                    [-]
% rsoi          [1x1]       radius Sphere of Influence      [km]
% mu            [1x1]       const                           [km^3/s^2]
%
% Output arguments:
% -----------------------------------------------------------------
% deltaT        [1x1]       delta time                      [s]
%
% CONTRIBUTORS:
%   Pier Francesco A. Bachini
%   Stefano Belletti
%   Chiara Giardini
%   Carolina Gómez Sánchez
%
% VERSION:
%   2024-01-10 latest

t0 = 0;

aBar = abs(a);
n = sqrt(mu/aBar^3);

theta = acos(((a*(1-e^2))/rsoi-1)/e);
F = 2*atanh(sqrt((e-1)/(e+1))*tan(theta/2));

t1 = (1/n) * (e*sinh(F) - F);

deltaT = t1 - t0;

end

