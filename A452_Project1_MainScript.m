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

% convert COEs to ECI
[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);
h_target_vector = cross(rECI.target,vECI.target);

%% CHASER initial conditions

% Get DCM Matrix from Target
QXx = QXx_from_rv_ECI(rECI.target,vECI.target);

% Get r.target in LVLH
rLVLH.target = QXx*rECI.target;

% Find r.chaser in ECI
rhoLVLH0 = [0;100;0]; % km
rLVLH.chaser = rLVLH.target - rhoLVLH0;
rECI.chaser = QXx' * rLVLH.chaser;

% Find velocity of Chaser in ECI
vECInorm.chaser = sqrt(mu/norm(rECI.chaser));
vECIdirect.chaser = cross(h_target_vector,rECI.chaser)/norm( cross(h_target_vector,rECI.chaser));
vECI.chaser = vECInorm.chaser * vECIdirect.chaser;

% Find LVLH velocity of Chaser and Target
vLVLH.chaser = QXx*vECI.chaser;
vLVLH.target = QXx*vECI.target;

% r and v from TLEs BOTH at mission time t0
[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);

% relative motion (shooting for 100km rho apart (LVLH))
[r_relx, v_relx, a_relx] = rva_relative(rECI.chaser,vECI.chaser,rECI.target,vECI.target); 

rho_rel = norm(r_relx); 
% disp(rho_rel + " km") % yay

%% Plot mission time T0 orbit(s)

% % % % Plot in ECI
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
p2 = plot3(newstate.target(end,1),newstate.target(end,2),newstate.target(end,3),'b','LineWidth',5);
p2.Marker = '*';

% Show CHASER at mission time t0
p3 = plot3(rECI.chaser(1),rECI.chaser(2),rECI.chaser(3),'k','LineWidth',5);
p3.Marker = '*';

% % % % Plot initial setup in LVLH

vbar0.chaser = rLVLH.target(2) - rLVLH.chaser(2);
rbar0.chaser = rLVLH.target(1) - rLVLH.chaser(1);

figure()
plot(0,0,'square','Linewidth',2) % target, center of LVLH frame
hold on
plot(vbar0.chaser,rbar0.chaser,'x','Linewidth',2) % chaser, in front of target in LVLH frame
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Relative Distance at Mission Start','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Chaser', 'interpreter','latex','Location', 'best')

%% First maneuver: Hop

% Two impulse
% burn on
% velocity initial and r_initial

% Parameters for trajectory burn
period = T.target; 
t = 3600*6; % for a 2-hour trajectory

drf = [0;40;0];
dr = r_relx;
dv0 = v_relx;
[hop1.dv0_PLUS_start_burn,hop1.dvf_MINUS_off_burn,hop1.deltaV,hop1.deltaV_afterBurn] = cw_twoimpulse(dr,drf,dv0, period,t);

disp("DeltaVs of Hop 1") 

% dv0_PLUS_start_burn % just display it
disp("Delta V: ")
disp(hop1.deltaV)
disp("Delta V After Burn: ")
disp(hop1.deltaV_afterBurn)

% Find new velocity vector of chaser to get ON the burn IN ECI
% current velocity in ECI == vECI.chaser
% hop1.dv0_chaser = -(vLVLH.chaser - dv0_PLUS_start_burn);
% hop1.vECIchaser0 = QXx*hop1.dv0_chaser;
% 
% 
% 
% 
breakpoint = 1; % WORKING HERE.




% Coast function over trajectory timespan
tspan = [0 t]; 
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state0_traj1.chaser = [rECI.chaser; hop1.vECIchaser0]; % situation at instant of Burn 1
state0_traj1.target = [rECI.target; vECI.target]; % situation at instant of Burn 1

% % % % % %
% Propogate TRAJECTORY 1 % WE NEED TO USE THE OTHER EOMMMMMMMMMS ODE CODE
% 12 inputs



% % % % % % time state h mu
%[hop1.time_chas, hop1.chaser] = ode45(@linearizedEOMs_std,tspan,state,options,h,mu);


hop1.r_chaser_after = hop1.chaser(:,1:3);
hop1.r_target_after = hop1.target(:,1:3);
hop1.v_chaser_after = hop1.chaser(:,4:6);
hop1.v_target_after = hop1.target(:,4:6);


% for i = 1:length(hop1.time_targ)
% QXx_loop1 = QXx_from_rv_ECI(hop1.r_target_after(i,1:3)',hop1.v_target_after(i,1:3)');
% rLVLH_afterhop1.target(i,1:3) = QXx_loop1 * hop1.r_target_after(i,1:3)';
% rLVLH_afterhop1.chaser(i,1:3) = QXx_loop1 * hop1.r_chaser_after(i,1:3)';
% end

% rLVLH_afterhop1.chaser = hop1.QXx * rECI_afterhop1.chaser';
% rECI_afterhop1.target = rECI_afterhop1.target';
% rECI_afterhop1.chaser = rECI_afterhop1.chaser';
% rLVLH_afterhop1.target = rLVLH_afterhop1.target';
% rLVLH_afterhop1.chaser = rLVLH_afterhop1.chaser';


figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(hop1.r_target_after(end,1),hop1.r_target_after(end,2),hop1.r_target_after(end,3),'r','LineWidth',2);
p2 = plot3(hop1.r_chaser_after(end,1),hop1.r_chaser_after(end,2),hop1.r_chaser_after(end,3),'k','LineWidth',2);
p1.Marker = 'square';
p2.Marker = '*';

% p2 = plot3(newstate.target(end,1),newstate.target(end,2),newstate.target(end,3),'b','LineWidth',5);
% p2.Marker = '*';
% 
% % Show CHASER at mission time t0
% p3 = plot3(rECI.chaser(1),rECI.chaser(2),rECI.chaser(3),'k','LineWidth',5);
% p3.Marker = '*';




% Extract r_bar values





% % % % % %
% Plot TRAJECTORY 1
% % % % % %
% figure()
% plot(0,0,'square','Linewidth',2) % target, center of LVLH frame
% hold on
% plot(time1.trajectory1,chaserhop1.traj_r,'x','Linewidth',2) % chaser, in front of target in LVLH frame
% xline(0)
% yline(0)
% 
% % Graph pretty 
% ylim padded 
% xlim padded 
% xLab = xlabel('Downrange [km]','Interpreter','latex'); 
% yLab = ylabel('Altitude [km]','Interpreter','latex'); 
% plotTitle = title('LVLH frame: Hop to 40km','interpreter','latex'); 
% set(plotTitle,'FontSize',14,'FontWeight','bold') 
% set(gca,'FontName','Palatino Linotype') 
% set([xLab, yLab],'FontName','Palatino Linotype') 
% set(gca,'FontSize', 9) 
% set([xLab, yLab],'FontSize', 14) 
% grid on 
% legend('Target','Chaser', 'interpreter','latex','Location', 'best')


% PLOT
% coast function 


% perform burn off
% velocity final 



