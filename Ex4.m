%% ex 1

close all
clear, clc

addpath("Functions\");

r1 = [-21800 37900 0];
r2 = [27300 27700 0];
deltaT = 15*3600 + 6*60 + 40;

mu = astroConstants(13);

orbitType = 0;

[A,P,E,ERROR,VI,VF,TPAR,THETA] = lambertMR(r1,r2,deltaT,mu,orbitType,0);

% propagation
[a, e, i, OM, om, theta] = rv2parorb(r1, VI, mu);
[theta_vect] = calculateThetaVect(mu, a, e, 200);
[rr, vv] = parorb2rv(a, e, i, OM, om, theta_vect, mu);

% plot
earthPlot;
plot3(rr(:,1), rr(:,2), rr(:,3))
axis equal, grid on
