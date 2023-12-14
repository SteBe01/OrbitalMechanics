function [plane_vect] = plane3points(v1, v2)

prod = cross(v1, v2);
plane_vect = prod ./ norm(prod);

% always up
if plane_vect(3) < 0
    plane_vect = -plane_vect;
end

end

