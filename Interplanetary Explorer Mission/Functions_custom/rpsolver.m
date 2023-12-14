function [rp] = rpsolver(v1, v2, mu)

a_minus = -mu/norm(v1)^2;
a_plus = -mu/norm(v2)^2;

delta = atan2((norm(cross(v1,v2)))/(norm(v1)*norm(v2)), (dot(v1,v2))/(norm(v1)*norm(v2)));

impact_param_minus = -a_minus/tan(delta/2);
impact_param_plus = -a_plus/tan(delta/2);

fun = @(rp) asin(1./(1+(rp*(norm(v1)^2))./mu))+asin(1./(1+(rp*(norm(v2)^2))./mu)) - delta;
rp = fzero(@(rp) fun(rp), (impact_param_plus+impact_param_minus)/2);

end

