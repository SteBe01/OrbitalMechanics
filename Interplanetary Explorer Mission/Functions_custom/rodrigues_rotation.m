function [v_rot] = rodrigues_rotation(v, delta, u)

% rotates a vector v of an angle delta around unit vector u
% (counter-clockwise)
%
% Usage
% [v_rot] = rodrigues_rotation(v,delta,u)
%
% Input arguments:
% ----------------------------------------------------------------
% v                 [3x1]   Vector                          [-]
% delta             [1x1]   Turning Angle                   [rad]
% u                 [3x1]   Unit Vector                     [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% v_rot             [3x1]   Rotated vector                  [km/s]

v_rot=v*cos(delta)+cross(u,v)*sin(delta)+u*dot(u,v)*(1-cos(delta));

end

