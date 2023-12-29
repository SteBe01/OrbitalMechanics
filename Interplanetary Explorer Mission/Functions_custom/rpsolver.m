function [rp, flag] = rpsolver(v1, v2, planetId)

% Perigee radius finder of hyperbola
%
% Usage
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

fun = @(rp) (delta - delta_fun(rp, v1, v2, mu)) * (rp > 0);

rp = fzero(@(rp) fun(rp), (impact_param_plus+impact_param_minus)/2);

if rp < astroConstants(planetId + 20)
    flag = 1;               % inside planet
else
    flag = 0;
end

end


function delta = delta_fun(rp, v_min, v_plus, mu)
    e_min = 1 + rp * norm(v_min)^2 / mu;
    e_plus = 1 + rp * norm(v_plus)^2 / mu;
    delta = asin(1 / e_min) + asin(1 / e_plus);
end

