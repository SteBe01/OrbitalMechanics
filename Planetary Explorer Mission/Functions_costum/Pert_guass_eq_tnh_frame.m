function [a_tnh] = Pert_guass_eq_tnh_frame(handleAcc,t,kep,mu,accArgs{:});
%Gauss' Equations in {t,n,h} reference frame system
at=a_tnh(1);
an=a_tnh(2);
ah=a_tnh(3);

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

da_dt=((2*v*a^2)/(mu))*at;
de_dt=(1/v)*(2*(e+cos(f))*at-(r/a)*sin(f)*an);
di_dt=(r*cos(f+om))*ah/h;
dOm_dt=(r*sin(f+om))*ah/(h*sin(i));
dom_dt=(1/(e*v))*(2*sin(f)*at+(2*e+(r/a)*cos(f))*an)-((r*sin(f+om)*cos(i))/(h*sin(i)))*ah;
df_dt=(h/r^2)+(1/(e*v))*(2*sin(f)*at+(2*e+(r/a)*cos(f))*an);
dM_dt=n-(b/(e*a*v))*(2*(1+((e^(2*r))/p))*sin(f)*at+(r/a)*cos(f)*an);
end