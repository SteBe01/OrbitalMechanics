function dy = eq_motion(t, s, acc_pert_fun, mu)
% Calculates the rate of change of keplerian elements on the t,n,h frame
%
% Usage:
% dy = eq_motion(t, s, acc_pert_fun, mu)
%
% Input arguments:
% ----------------------------------------------------------------
% t             [Nx1]   time       [-]
% s             [1x6]   Keplerian Elements Vector   [-]
% acc_pert_fun  [1x3]   Perturbing Accelerations    [km/s^2]
% mu            [1x1]   Planetary constant          [km^3/s^2]
%
% Output arguments:
% -----------------------------------------------------------------
% dy            [6x1] Rate of change of keplerian elements    [km/s^2]

    a = s(1);
    e = s(2);
    i = s(3);
    OM = s(4);
    om = s(5);
    th = s(6);

    acc_pert_vec = acc_pert_fun(t, s);
    acc_t = acc_pert_vec(1); acc_n = acc_pert_vec(2); acc_h = acc_pert_vec(3);

    p = a*(1-e^2);
    r = p/(1+e*cos(th));
    v = sqrt(2*mu/r - mu/a);
    h = sqrt(p*mu);


    dy(1) = 2*a^2*v/mu * acc_t;
    dy(2) = 1/v * (2*(e+cos(th))*acc_t - r/a*sin(th)*acc_n);
    dy(3) = r*cos(th+om)/h * acc_h;
    dy(4) = r*sin(th+om)/(h*sin(i)) * acc_h;
    dy(5) = 1/(e*v) * (2*sin(th)*acc_t + (2*e + r/a*cos(th))*acc_n) - ...
        (r*sin(th+om)*cos(i))/(h*sin(i))*acc_h;
    dy(6) = h/r^2 - 1/(e*v) * (2*sin(th)*acc_t + (2*e + r/a*cos(th))*acc_n);

    dy = dy';

end