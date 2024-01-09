function [dv] = completeInterplanetaryMS(t1, t2, t3, code1, code2, code3)

% Function for interplanetary dv calculator - optimization only (multi
% start)
%
% Usage
% [dv] = completeInterplanetaryGS(t1, t2, t3, code1, code2, code3)
%
% Input arguments:
% ----------------------------------------------------------------
% t1            [1x1]       time instant for departure      [mjd200]
% t2            [1x1]       time instant for flyby          [mjd200]
% t3            [1x1]       time instant for arrival        [mjd200]
% code1         [1x1]       departure body Id               [-]
% code2         [1x1]       flyby body Id                   [-]
% code3         [1x1]       arrival body Id                 [-]
%
% -----------------------------------------------------------------
% Output arguments:
% 
% dv            [1x1]       total dv                        [km/s]

if t1 < date2mjd2000([2028 01 01 0 0 0]) || t3 > date2mjd2000([2058 01 01 0 0 0])
    dv = 1e7;
    return
end

tof1 = (t2 - t1) * 24 * 60 * 60;
tof2 = (t3 - t2) * 24 * 60 * 60;

if tof1 <= 0 || tof2 <= 0
    dv = 1e7;
    return
end

if code1 <= 10
    [departure.kep, ksun] = uplanet(t1, code1);
else
    [departure.kep, ~] = ephNEO(t1, code1);  
end

if code2 <= 10
    [flyby.kep, ~] = uplanet(t2, code2);
else
    [flyby.kep, ~] = ephNEO(t2, code2);  
end
if code3 <= 10
    [arrival.kep, ~] = uplanet(t3, code3);
else
    [arrival.kep, ~] = ephNEO(t3, code3);  
end

[departure.r0, departure.v0] = kep2car([departure.kep ksun]);
[flyby.r0, flyby.v0] = kep2car([flyby.kep ksun]);
[arrival.r0, arrival.v0] = kep2car([arrival.kep ksun]);

orbitType = 0;
nOrbits = 0;
[A_1, P_1, E_1, ERROR_1, VI_1, VF_1, TPAR_1, THETA_1] = lambertMR(departure.r0, flyby.r0, tof1, ksun, orbitType, nOrbits);
[A_2, P_2, E_2, ERROR_2, VI_2, VF_2, TPAR_2, THETA_2] = lambertMR(flyby.r0, arrival.r0, tof2, ksun, orbitType, nOrbits);

flyby.v_inf_minus = VF_1 - flyby.v0;
flyby.v_inf_plus = VI_2 - flyby.v0;

if norm(flyby.v_inf_minus) > 1e8 || norm(flyby.v_inf_plus) > 1e8
    dv = 1e7;
    return
end

if code2 + 10 > 19
    dv = 1e7;
    return
end
flyby.mu = astroConstants(code2 + 10);
% rp = 1e3;
rp = rpsolver(flyby.v_inf_minus, flyby.v_inf_plus, code2); % posizione giusta?
if rp < astroConstants(code2 + 20)
    dv = 1e7;
    return
end

if rp > 48.2e6
    dv = 1e7;
    return
end

dv1 = norm(VI_1 - departure.v0);
dv2 = abs(sqrt((2*flyby.mu/rp)+norm(flyby.v_inf_plus)^2)-sqrt((2*flyby.mu/rp)+norm(flyby.v_inf_minus)^2));
dv3 = norm(arrival.v0 - VF_2);
dv = dv1 + dv2 + dv3;


end

