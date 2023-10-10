%{
 AERO 452: Spaceflight Dynamics II
Group Project #1
Group #21
Authors: Suren Sanai and Justin Self 
Main Script here. 
All functions are separate .m files and MUST be stored in the GitHub local
folder.
%}

% Housekeeping
clear all; close all; clc;

% Constants
mu = 398600; % km3/s2
r_earth = 6378; % km

%% Target Information

% COEs from Heavens-Above
%{

INTELSAT 32E - Orbit
Sat ID: 41945

Epoch (UTC):	04 October 2023 20:36:45
Eccentricity:	0.0000178
inclination:	0.0182°
perigee height:	35786 km
apogee height:	35787 km
right ascension of ascending node:	167.1380°
argument of perigee:	309.8738°
revolutions per day:	1.00270360
mean anomaly at epoch:	162.2886°
orbit number at epoch:	2434
%}

% NEXT STEPS
% COEs to r,v (Justin)
% Put our chaser vehicle 100km away (rho, LVLH) from Target
% use universal var to do this. 


%% TARGET initial conditions

RAAN.target = 167.1380; % deg
inc.target = 0.0182; % deg
w.target = 309.8738; % deg
alt.target = (35786 + 35787) / 2; % altitude
r.target = alt.target - r_earth; % km
ecc.target = 0.0000178;
theta.target = 162.2886; % deg

h.target = findh(r.target,mu,ecc.target,theta.target);

[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);


%% CHASER initial conditions

RAAN.chaser = 167.1380; % deg
inc.chaser = 0.0182; % deg
w.chaser = 309.8738; % deg
alt.chaser = (35786 + 35787) / 2; % altitude
r.chaser = alt.target - r_earth; % km
ecc.chaser = 0.0000178;
theta.chaser = theta.target-.2; % deg

h.chaser = findh(r.chaser,mu,ecc.chaser,theta.chaser);

[rECI.chaser,vECI.chaser] = r_and_v_from_COEs(RAAN.chaser,inc.chaser,w.chaser,h.chaser,ecc.chaser,theta.chaser);

% relative motion
[r_relx, v_relx, a_relx] = rva_relative(rECI.chaser,vECI.chaser,rECI.target,vECI.target); 

rho_rel = norm(r_relx); 

disp(rho_rel + " km")

