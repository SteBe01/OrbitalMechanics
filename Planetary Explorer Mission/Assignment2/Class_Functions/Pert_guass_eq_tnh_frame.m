function [acc_tnh] = Pert_guass_eq_tnh_frame(t,kep,acc_pert_fun,mu)
%Gauss' Equations in {t,n,h} reference frame system

%Acceleration elements
acc_pert_vec = acc_pert_fun(t, kep);
acc_t = acc_pert_vec(1); 
acc_n = acc_pert_vec(2); 
acc_h = acc_pert_vec(3);

%keplerian elements
a=kep(1); %[km]
e=kep(2);
i=kep(3); %[rad]
Om=kep(4); %[rad]
om=kep(5); %[rad]
f=kep(6); %[rad]

%orbit values 
b=a*sqrt(1-e^2); %[km]
p=(b^2)/a; %[km]
r=p/(1+e*cos(f)); %[km]
v=sqrt((2*mu/r)-(mu/a)); %[km/s]
n=sqrt(mu/(a^3));
h=n*a*b;

%guass planetary equations
da_dt=((2*v*a^2)/(mu))*acc_t;
de_dt=(1/v)*(2*(e+cos(f))*acc_t-(r/a)*sin(f)*acc_n);
di_dt=(r*cos(f+om))*acc_h/h;
dOm_dt=(r*sin(f+om))*acc_h/(h*sin(i));
dom_dt=(1/(e*v))*(2*sin(f)*acc_t+(2*e+(r/a)*cos(f))*acc_n)-((r*sin(f+om)*cos(i))/(h*sin(i)))*acc_h;
df_dt=(h/r^2)+(1/(e*v))*(2*sin(f)*acc_t+(2*e+(r/a)*cos(f))*acc_n);
dM_dt=n-(b/(e*a*v))*(2*(1+((e^(2*r))/p))*sin(f)*acc_t+(r/a)*cos(f)*acc_n);

acc_tnh=[da_dt  de_dt  di_dt  dOm_dt  dom_dt  df_dt]';

end