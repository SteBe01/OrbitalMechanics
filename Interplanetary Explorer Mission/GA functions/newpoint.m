function [v_new] = newpoint(v_impact, ref_point, impact_param, hyp_planeVect)

hyp_planeVect = hyp_planeVect ./ norm(hyp_planeVect);

v_impact = v_impact - ref_point;

v_impact_module = norm(v_impact);
alpha = asin(impact_param/v_impact_module);

v_new_versor = v_impact*cos(alpha) + cross(hyp_planeVect, v_impact)*sin(alpha) + hyp_planeVect*dot(hyp_planeVect, v_impact)*(1 - cos(alpha));
v_new_versor = v_new_versor./norm(v_new_versor);

v_new = sqrt(v_impact_module^2 - impact_param^2) .* v_new_versor;

v_new = v_new + ref_point;

end

