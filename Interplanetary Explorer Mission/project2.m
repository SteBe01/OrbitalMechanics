

addpath("Functions\")

clear, clc
close all

% data - goup 2346
%Departure: Saturn
%Flyby: Jupiter
%Arrival: Asteriod N.79
%Earliest Departure: 00:00:00 01/01/2028
%Latest Arrival: 00:00:00 01/01/2058

%in heliocentric region
rsoirsh=0.0382; %for Saturn
rsoirjh=0.0620; %for Jupiter
rsoirah=0.0620; %for Asteriod N.79 using Jupiter's Value

%in planetocentric region
rsoirSp=906.9; %for Saturn
rsoirJp=675.1; %for Jupiter
rsoirAp=675.1; %for Asteriod N.79 using Jupiter's Value

%Values defined by us for flyby
%radius of percentre

%Excess velocity
