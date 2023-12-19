function [a_rsw] = Pert_guass_eq_rsw_frame(handleAcc,t,kep,mu,accArgs{:});
%Gauss' Equations in {r,s,w} reference frame system
ar=a_rsw(1);
as=a_rsw(2);
aw=a_rsw(3);

a=kep(1); %[km]
e=kep(2);
i=kep(3); %[rad]
Om=kep(4); %[rad]
om=kep(5); %[rad]
f=kep(6); %[rad]

b=a*sqrt(1-e^2); %[km]
p=(b^2)/a; %[km]
r=p/(1+e*cos(f)); %[km]
v=sqrt((2*mu/r)-(mu/a)); %[km/s]
n=sqrt(mu/(a^3));
h=n*a*b;

da_dt=(2*a^2)*(e*sin(f)*ar+(p/r)*as)/h;
de_dt=(1/h)*(p*sin(f)*ar+((p+r)*cos(f)+r*e)*as);
di_dt=(r*cos(f+om))*aw/(h);
dOm_dt=(r*sin(f+om))*aw/(h*sin(i));
dom_dt=(1/(h*e))*(-p*cos(f)*ar+(p+r)*sin(f)*as)-(r*sin(f+om)*cos(i))*aw/(h*sin(i));
df_dt=(h/r^2)+(1/(h*e))*(p*cos(f)*ar-(p+r)*sin(f)*as);
dh_dt=r*as;
dM_dt=n+(b/(a*h*e))*((p*cos(f)-2*r*e)*ar-(p+r)*sin(f)*as);
end