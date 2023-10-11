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
T.target = 2*pi*sqrt(r.target^3/mu);

h.target = findh(r.target,mu,ecc.target,theta.target);


[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);
h_target_vector = cross(rECI.target,vECI.target);
%% CHASER initial conditions

RAAN.chaser = 167.1380; % deg
inc.chaser = 0.0182; % deg
w.chaser = 309.8738; % deg
alt.chaser = (35786 + 35787) / 2; % altitude
r.chaser = alt.target - r_earth; % km
ecc.chaser = 0.0000178;
theta.chaser = theta.target -.1948226; % deg
T.chaser = 2*pi*sqrt(r.chaser^3/mu);
h.chaser = findh(r.chaser,mu,ecc.chaser,theta.chaser);


% Get DCM Matrix from Target
QXx = QXx_from_rv_ECI(rECI.target,vECI.target);

% Get r.target in LVLH
rLVLH.target = QXx*rECI.target;

% Find r.chaser in ECI
rhoLVLH0 = [0;100;0]; % km
rLVLH.chaser = rLVLH.target - rhoLVLH0;
rECI.chaser = QXx' * rLVLH.chaser;

% Find velocity of Chaser in ECI
vECIchaser_norm = sqrt(mu/norm(rECI.chaser));
vECI_chaser_direct = cross(h_target_vector,rECI.chaser)/norm( cross(h_target_vector,rECI.chaser));
vECI.chaser = vECIchaser_norm * vECI_chaser_direct;


% r and v from TLEs BOTH at mission time t0

[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);

% relative motion (shooting for 100km rho apart (LVLH))
[r_relx, v_relx, a_relx] = rva_relative(rECI.chaser,vECI.chaser,rECI.target,vECI.target); 

rho_rel = norm(r_relx); 

disp(rho_rel + " km") % yay

%% Plot mission time T0 orbit(s)
tspan = [0 T.target]; % seconds = 1 day
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state0.target = [rECI.target,vECI.target];

% Propogate 1 period TARGET
[newtime.target, newstate.target] = ode45(@coast_ODE,tspan,state0.target,options,mu);

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(newstate.target(:,1),newstate.target(:,2),newstate.target(:,3),'r','LineWidth',2);
p2 = plot3(newstate.target(end,1),newstate.target(end,2),newstate.target(end,3),'k','LineWidth',5);
p2.Marker = '*';

% Show CHASER at mission time t0
plot3(rECI.chaser(1),rECI.chaser(2),rECI.chaser(3),'*','LineWidth',5)


