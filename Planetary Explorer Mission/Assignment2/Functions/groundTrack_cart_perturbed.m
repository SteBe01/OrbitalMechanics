function [alpha, delta, lon, lat] = groundTrack_cart_perturbed(y0, t_vect, mu_E, theta0_g, om_E, Re, J2, A_M, cD)

% alpha, delta, lon and lat calculator for perturbed orbits
%
% usage:
% [alpha, delta, lon, lat] = groundTrack_cart_perturbed(y0, t_vect, mu_E, theta0_g, om_E, Re, J2, A_M, cD)
%
% Input arguments:
% ----------------------------------------------------------------
% y0            [6x1]       y0 = [r0'; v0']                 [-]
% t_vect        [1xN]       time vector                     [s]
% mu_E          [1x1]       mu Earth                        [km^3/s^2]
% theta0_g      [1x1]       theta0 of Earth                 [deg]
% om_E          [1x1]       angular velocity of Earth       [deg/h]
% TO BE FINISHED
% 
% Output arguments:
% -----------------------------------------------------------------
% alpha         [Nx1]       alpha                           [deg]
% delta         [Nx1]       delta                           [deg]
% lon           [Nx1]       lon                             [deg]
% lat           [Nx1]       lat                             [deg]
% 
% CONTRIBUITORS:
% Pier Francesco A. Bachini
% Stefano Belleti
% Chiara Giardini
% Carolina Gómez Sánchez

options = odeset('RelTol', 1e-13, 'AbsTol', 1e-14);
[~, Y] = ode113(@(t,y) ode_2bp_perturbed( t, y, mu_E, Re, J2, deg2rad(om_E) / 3600, A_M, cD), t_vect, y0, options);

r = vecnorm(Y(:,1:3),2,2);

alpha = atan2(Y(:,2), Y(:,1)) * (180/pi);
delta = asin(Y(:,3)./r) * (180/pi);
om_E = om_E / 3600;                                                         % from h to seconds
lon = alpha - (ones(length(t_vect),1)*theta0_g + om_E*t_vect');
lat = delta;

lon = lon + ones(length(lon),1)*180;

for i=1:length(t_vect)
    if lon(i)<0
        while lon(i) < -360
            lon(i) = lon(i) + 360;
        end
    end
    if lon(i)<0
        lon(i) = lon(i) + 360;
    end
end
lon = lon - ones(length(lon),1)*180;

end
