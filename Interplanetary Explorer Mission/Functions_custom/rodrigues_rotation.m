function [v_rot] = rodrigues_rotation(v,delta,u)

% rotate a vector 𝐯 an angle 𝛿 around unit vector 𝐮 (counter-clockwise):
%
% Usage
% [v_rot] = rodrigues_rotation(v,delta,u)
%
% Input arguments:
% ----------------------------------------------------------------
% v                 [3x1]   Velocity Vector                            [km/s]
% delta             [1x1]   Turning Angle                              [rad]
% u                 [3x1]   Unit Vector normal to Perifocal frame      [-]
% 
% Output arguments:
% -----------------------------------------------------------------
% a                 [3x1]   Velocity Vector                            [km/s]

v_rot=v*cos(delta)+cross(u,v)*sin(delta)+u*dot(u,v)*(1-cos(delta));

end