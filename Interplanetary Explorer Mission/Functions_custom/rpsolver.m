function [rp, flag] = rpsolver(v1, v2, planetId)

% Perigee radius finder of hyperbola
%
% [rp, flag] = rpsolver(v1, v2, planetId)
%
% Input arguments:
% ----------------------------------------------------------------
% v1            [1x1]       velocity -inf           [km/s]
% v2            [1x1]       velocity +inf           [km/s]
% plnetId       [1x1]       planet Id               [-]
%
% -----------------------------------------------------------------
% Output arguments:
% 
% rp            [1x1]       perigee radius          [km]
% flag          [1x1]       exit flag               [km]


mu = astroConstants(planetId + 10);

a_minus = -mu/norm(v1)^2;
a_plus = -mu/norm(v2)^2;

% delta = atan2((norm(cross(v1,v2)))/(norm(v1)*norm(v2)), (dot(v1,v2))/(norm(v1)*norm(v2)));
delta = real(acos(dot(v1, v2)/(norm(v1)*norm(v2))));

impact_param_minus = -a_minus/tan(delta/2);
impact_param_plus = -a_plus/tan(delta/2);

fun = @(rp) (asin(1./(1+(rp*(norm(v1)^2))./mu))+asin(1./(1+(rp*(norm(v2)^2))./mu)) - delta) * rp > 0;
rp = fzero(@(rp) fun(rp), impact_param_plus+impact_param_minus)/2;

if rp < astroConstants(planetId + 20)
    flag = 0;               % inside planet
else
    flag = 1;
end

end

