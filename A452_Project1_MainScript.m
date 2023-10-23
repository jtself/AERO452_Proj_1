%{
 AERO 452: Spaceflight Dynamics II
Group Project #1
Group #21
Authors: Suren Sanai and Justin Self 
Main Script.
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
% COEs to r,v
% Put our chaser vehicle 100km away (rho, LVLH) from Target
% use universal var to do this. 

%% TARGET initial conditions

RAAN.target = 167.1380; % deg
inc.target = 0.0182; % deg
w.target = 309.8738; % deg
alt.target = (35786 + 35787) / 2; % altitude
r.target = alt.target + r_earth; % km
ecc.target = 0.0000178;
theta.target = 162.2886; % deg
T.target = 2*pi*sqrt(r.target^3/mu);

h.target = findh(r.target,mu,ecc.target,theta.target);

[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);
h_target_vector = cross(rECI.target,vECI.target);
%% CHASER initial conditions
% Written by JS and SS
% Approved by JS

% Chaser must be 100km away EXACTLY from the target at mission t0. 

% Get DCM Matrix from Target
QXx = QXx_from_rv_ECI(rECI.target,vECI.target);

% Get r.target in LVLH
rLVLH.target = QXx*rECI.target;

% Find r.chaser in ECI
rhoLVLH0 = [0;100;0]; % km FORCING FUNCTION per mission parameter
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

% relative motion
[r_relx, v_relx, a_relx] = rva_relative(rECI.chaser,vECI.chaser,rECI.target,vECI.target); 

rho_rel = norm(r_relx); 

% Print check
disp("* CHECK *  Relative distance between target and chaser at mission t0 is: " + rho_rel + " km")

%% Plot mission time T0 orbit(s)
% Written by JS
% Approved by SS

tspan = [0 T.target]; % seconds = 1 day
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state0.target = [rECI.target,vECI.target];

% Propogate 1 period TARGET
[newtime0.target, newstate0.target] = ode45(@coast_ODE,tspan,state0.target,options,mu);

figure()
   h1 = gca;
   earth_sphere(h1)
   hold on

% % TARGET at mission start time, t0

% Plot orbit path
p1 = plot3(newstate0.target(:,1),newstate0.target(:,2),newstate0.target(:,3),'r','LineWidth',1); 

% Plot Target position in orbit
p2 = plot3(newstate0.target(end,1),newstate0.target(end,2),newstate0.target(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER at mission time t0
pc = plot3(rECI.chaser(1),rECI.chaser(2),rECI.chaser(3),'x','LineWidth',2);
pc.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('Spacecrafts A and B at Mission $t_0$','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('','Common orbit','Target (A)','Chaser (B)', 'interpreter','latex','Location', 'best')

%% Plot mission t0 configuration in LVLH
% Written by JS and SS
% Approved by JS 10/21/23

vbar0.chaser = rLVLH.target(2) - rLVLH.chaser(2);
rbar0.chaser = rLVLH.target(1) - rLVLH.chaser(1);

figure()
pt = plot(0,0,'square','Linewidth',2); % target, center of LVLH frame
pt.Color = 'blue';
hold on
pc = plot(vbar0.chaser,rbar0.chaser,'x','Linewidth',2);
pc.Color = 'k';
% chaser, in front of target in LVLH frame
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Relative Distance at Start of Mission ($t_0$)','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target (A)','Chaser (B)', 'interpreter','latex','Location', 'best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FIRST RENDEZVOUS MANEUVER: 100km to 40km hop
% Written by JS and SS
% Approved by JS 10/21/23

period = T.target; 

% Choose trajectory travel time.
hours = 12; % This minimized delta-v
t = 3600*hours;

% Where you want to end up
drf = [0;40;0];

% What the current relative distance and velocity (LVLH) is
dr = r_relx;
dv0 = v_relx;

% Call function to find instantaneous dv burn (start of trajectory)
[hop1.dv0_plus,hop1.DV_0,hop1.dV_F,hop1.DV_total] = cw_twoimpulse(dr,drf,dv0,period,t);

%{
OUTPUTS: 
dv0_plus = chaser initial burn dv to get onto transfer trajectory. 
DV_0 = total Delta-v used up in first burn
dv_f = second burn to get off transfer trajectory
DV_total = total DV in full manuever (both burns)
%}

%{
Use linearized EOMs to hop. 
state = [delx
        dely
        delz
        delxdot
        delydot
        delzdot]
%}

% % Propogate transfer trajectory

tspan = [0 t]; % length of trajectory flight
dv = hop1.dv0_plus; % This is the new relative velocity start burn
state = [dr;dv;rECI.target;vECI.target];
[timenew,statenew] = ode45(@linearizedEOMs_std,tspan,state,options,h.target,mu);

% Extract data after ODE
hop1.rECI_target_data = [statenew(:,7),statenew(:,8), statenew(:,9)];
hop1.vECI_target_data = [statenew(:,10),statenew(:,11), statenew(:,12)];
relativePosition_Vect = [statenew(:,1),statenew(:,2),statenew(:,3)]; % since z is zero whole time
relativeVelocity_Vect = [statenew(:,4),statenew(:,5),statenew(:,6)]; 

relativePosition = relativePosition_Vect(end,1:3);
hop1.verifyDist = norm(relativePosition); % should be 40km

disp("* CHECK *  Relative distance between target and chaser after maneuver is: " + hop1.verifyDist + " km")

%% Plot trajectory: 100km to 40km hop
% Written by JS and SS
% Approved by JS 10/21/23

figure()
% target, center of LVLH frame
pc = plot(0,0,'square','Linewidth',2);
pc.Color = 'b';
hold on
% Hop trajectory
plot(relativePosition_Vect(:,2),relativePosition_Vect(:,1),'r','LineWidth',2)

% Chaser position after hop
% Plot
p1 = plot(relativePosition_Vect(end,2),relativePosition_Vect(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 100 km to 40 km hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% NEXT: Plot this whole (FIRST) maneuver in ECI
% Written by JS and SS
% Approved by JS 10/21/23

% Convert LVLH state data of the chaser on the first hop to ECI

for i = 1:length(timenew)
    hop1QXx = QXx_from_rv_ECI(hop1.rECI_target_data(i,:)',hop1.vECI_target_data(i,:)');
    hop1.rECI_chaser_data(i,:) = (hop1QXx' * relativePosition_Vect(i,:)') + hop1.rECI_target_data(i,:)';
    hop1.vECI_chaser_data(i,:) = (hop1QXx' * relativeVelocity_Vect(i,:)') + hop1.vECI_target_data(i,:)';
end

figure()
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
% Orbital path
p1 = plot3(hop1.rECI_target_data(:,1),hop1.rECI_target_data(:,2),hop1.rECI_target_data(:,3),'--k','LineWidth',1);
% Target location in orbit
p2 = plot3(hop1.rECI_target_data(end,1),hop1.rECI_target_data(end,2),hop1.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER at mission time t0
% Chaser orbit
p3 = plot3(hop1.rECI_chaser_data(:,1),hop1.rECI_chaser_data(:,2),hop1.rECI_chaser_data(:,3),'r','LineWidth',1);
p3.Color = 'r';
% Chaser location at end of hop maneuver
p4 = plot3(hop1.rECI_chaser_data(end,1),hop1.rECI_chaser_data(end,2),hop1.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('100 km to 40 km hop maneuver in ECI','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('','Target orbit','Target','Chaser orbit during hop','Chaser', 'interpreter','latex','Location', 'best')

%% MANEUVER 1: 100 KM TO 40 KM HOP SUMMARY 
% Written by JS 10/21/23

disp(" ")
disp("--------- Hop: 100km to 40km ---------") 

missionDV.hop1 = hop1.DV_total*1000;
% Print DV_total from this maneuver
disp("Total Delta V for this maneuver is: " + hop1.DV_total*1000 + " m/s") % total mission DV to date

% Keep track of mission time
MET.hop1 = hours; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET (just after hop 1) = " + MET.hop1 + " hours")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NEXT: Burn off this hop trajectory onto v-bar in a static fashion FOOTBALL
% Written originally by SS
% Approved by JS 10/21/23

% NOTE: we do not "burn off" the hop trajectory. We simply jump right into
% the football orbit. NEED TO HOLD in the FOOTBALL ORBIT fo x hours.


n.Target = (2*pi)/T.target; % 1/seconds
b = relativePosition_Vect(end,2)/2; % relative position in the Y-DIRECTION after 100 - 40 km hop
Football1.xdot0  = b * n.Target; % altitude direction
dv_FootBall1_LVLH = [Football1.xdot0 0 0]; 

% Relative velocity after the trajectory hop burn == preBurn Football
Football1.dvChaser0 = relativeVelocity_Vect + dv_FootBall1_LVLH;

delta_v_intoFB1 = norm(Football1.dvChaser0) / 1000; % in m/s now

% disp(" ")
% disp("------------------ ")
% disp("Delta V to burn into Football Manuever 1: " + delta_v_intoFB1 + " m/s")

% Position after hop 1 / beginning of football
dr_t = relativePosition(end,1:3);

% Relative velocity 
dv_t = dv_FootBall1_LVLH;
dr0 = dr_t;
dv0 = dv_t;

% t = T.target/50000; % DONT NEED THIS

% Preallocate vectors
Football1.relativePosition = zeros(1,3);
Football1.relativeVelocity = zeros(1,3);

% ----- ALTERNATE METHOD USING CW MATRIX ---------
% for i = 1:50000
% [dr_t,dv_t] = CW_Matrix(t,dr0',dv0',n.Target);
% Football1.relativePosition(i,1:3) = dr_t;
% Football1.relativeVelocity(i,1:3) = dv_t;
% dr0 = dr_t';
% dv0 = dv_t';
% end

% This is old; the inputs to the funct were backwards % 
%{
% relative motion
[Football1.r_relx0, Football1.v_relx0, Football1.a_relx0] = rva_relative(...
    hop1.rECI_chaser_data(end,1:3)',hop1.vECI_chaser_data(end,1:3)',...
    hop1.rECI_target_data(end,1:3)',hop1.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse
%}

% Inputs for relative position calculation
rA = hop1.rECI_target_data(end,1:3)';
vA = hop1.vECI_target_data(end,1:3)';
rB = hop1.rECI_chaser_data(end,1:3)';
vB = hop1.vECI_chaser_data(end,1:3)';

% relative motion
[Football1.r_relx0, Football1.v_relx0, Football1.a_relx0] = rva_relative(rA,vA,rB,vB);

%% Hold on the football orbit for ONE TARGET ORBITAL PERIOD

% Length of hold
tspan = [0 T.target]; % length of trajectory flight == one target orbit

% Prepare to perform football trajectory insertion
Football1.dr0 = relativePosition;
Football1.dv0 = relativeVelocity_Vect(end,1:3) + dv_FootBall1_LVLH;
dr = relativePosition(end,1:3);
dv = Football1.dvChaser0;
state_FB1 = [dr0';dv0';...
    hop1.rECI_target_data(end,1:3)';hop1.vECI_target_data(end,1:3)'];
[Football1.timenew,Football1.statenew] = ode45(@linearizedEOMs_std,tspan,state_FB1,options,h.target,mu);


% Extract data after ODE
Football1.rECI_target_data = [Football1.statenew(:,7),Football1.statenew(:,8), Football1.statenew(:,9)];
Football1.vECI_target_data = [Football1.statenew(:,10),Football1.statenew(:,11), Football1.statenew(:,12)];
Football1.relativePosition_Vect = [Football1.statenew(:,1),Football1.statenew(:,2),Football1.statenew(:,3)]; % since z is zero whole time
Football1.relativeVelocity_Vect = [Football1.statenew(:,4),Football1.statenew(:,5),Football1.statenew(:,6)]; 

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on

% Plot Football trajectory
plot(Football1.relativePosition_Vect(:,2),Football1.relativePosition_Vect(:,1),'r','LineWidth',2)

% Chaser position after hop/at the football start
% Plot
p1 = plot(Football1.relativePosition_Vect(end,2),Football1.relativePosition_Vect(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('HOLD 1: 40 km - 20 km football orbit (LVLH frame)','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Chaser Trajectory (football orbit)', 'Chaser initial/final position','interpreter','latex','Location', 'best')


%% Plot both maneuvers in LVLH
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
xline(0)
yline(0)

% Hop 1 trajectory
plot(relativePosition_Vect(:,2),relativePosition_Vect(:,1),'--k','LineWidth',1)

% Plot Football trajectory
plot(Football1.relativePosition_Vect(:,2),Football1.relativePosition_Vect(:,1),'r','LineWidth',2)

% Chaser position after hop/at the football start
% Plot
p1 = plot(Football1.relativePosition_Vect(end,2),Football1.relativePosition_Vect(end,1),'x','LineWidth',2);
p1.Color = 'k';

% % CHECK football direction:  DONT NEED THIS ON FINAL VERSION; CHECK ONLY
% p1 = plot(Football1.relativePosition_Vect(100,2),Football1.relativePosition_Vect(100,1),'x','LineWidth',5);
% p1.Color = 'k';

% Graph pretty 
ylim([-50 50]) 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: First two maneuvers','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','','','100 km to 40 km hop','Football hold', 'Chaser after hold','interpreter','latex','Location', 'best')

%% NEXT: Plot the first FOOTBALL maneuver in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Football1.timenew)
    FB1QXx = QXx_from_rv_ECI(Football1.rECI_target_data(i,:)',Football1.vECI_target_data(i,:)');
    Football1.rECI_chaser_data(i,:) = (FB1QXx' * Football1.relativePosition_Vect(i,:)') + Football1.rECI_target_data(i,:)';
    Football1.vECI_chaser_data(i,:) = (FB1QXx' * Football1.relativeVelocity_Vect(i,:)') + Football1.vECI_target_data(i,:)';
end

figure()
   h1 = gca;
   earth_sphere(h1)
   hold on

% target orbit trail
p1 = plot3(Football1.rECI_target_data(:,1),Football1.rECI_target_data(:,2),Football1.rECI_target_data(:,3),'--k','LineWidth',1);
% Target location
p2 = plot3(Football1.rECI_target_data(end,1),Football1.rECI_target_data(end,2),Football1.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER after football hold
p3 = plot3(Football1.rECI_chaser_data(:,1),Football1.rECI_chaser_data(:,2),Football1.rECI_chaser_data(:,3),'r','LineWidth',1);
p4 = plot3(Football1.rECI_chaser_data(end,1),Football1.rECI_chaser_data(end,2),Football1.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('HOLD 1: Football orbit in ECI','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('','Target Orbit','Target','Chaser Path', 'Chaser', 'interpreter','latex','Location', 'best')

% <<<<<<< HEAD
%% SECOND HOP - from 40 km to 1 km
% =======
%% MANEUVER 2: FOOTBALL HOLD SUMMARY 
% Written by JS 10/21/23

disp(" ")
disp("--------- Hold 1: Football at 20 km - 40 km ---------") 
% >>>>>>> main

% Choose trajectory travel time.
hours = 6; % your choice
t = 3600*hours;

% Where you want to end up
Hop2.drf = [0;1;0];

% relative motion for two impulse inputs
[Hop2.r_relx0, Hop2.v_relx0, Hop2.a_relx0] = rva_relative(...
    Football1.rECI_chaser_data(end,1:3)',Football1.vECI_chaser_data(end,1:3)',...
    Football1.rECI_target_data(end,1:3)',Football1.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse


% What the current relative distance and velocity (LVLH) is
Hop2.dr = Football1.relativePosition(end,1:3)';
Hop2.dv0 = Football1.relativeVelocity(end,1:3)'; 
% Call function to find instantaneous dv burn (start of trajectory)
[Hop2.dv0_PLUS_start_burn,Hop2.dvf_MINUS_off_burn,Hop2.deltaV,Hop2.deltaV_afterBurn] = cw_twoimpulse(Hop2.dr,Hop2.drf,Hop2.dv0,period,t);

disp(" ")
disp("Hop 2: 40km to 1km") 

% dv0_PLUS_start_burn % just display it
disp("Hop 2 Delta V: ")
disp(Hop2.deltaV)
disp("Hop 2Delta V After Burn: ")
disp(Hop2.deltaV_afterBurn)



% SECOND HOP - from 40 km to 1 km

% Choose trajectory travel time.

t = T.target;

% Where you want to end up
Hop2.drf = [0;1;0];

% relative motion for two impulse inputs
[Hop2.r_relx0, Hop2.v_relx0, Hop2.a_relx0] = rva_relative(...
    Football1.rECI_chaser_data(end,1:3)',Football1.vECI_chaser_data(end,1:3)',...
    Football1.rECI_target_data(end,1:3)',Football1.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse


% What the current relative distance and velocity (LVLH) is
Hop2.dr = Football1.relativePosition(end,1:3)';
Hop2.dv0 = Football1.relativeVelocity(end,1:3)'; 
% Call function to find instantaneous dv burn (start of trajectory)
[Hop2.dv0_PLUS_start_burn,Hop2.dvf_MINUS_off_burn,Hop2.deltaV,Hop2.deltaV_afterBurn] = cw_twoimpulse(Hop2.dr,Hop2.drf,Hop2.dv0,period,t);

disp(" ")
disp("Hop 2: 40km to 1km") 

% dv0_PLUS_start_burn % just display it
disp("Delta V: ")
disp(Hop2.deltaV)
disp("Delta V After Burn: ")
disp(Hop2.deltaV_afterBurn)


%% Plot trajectory: 40km to 1km hop in LVLH

%{
Use linearized EOMs to hop. 
state = [delx
        dely
        delz
        delxdot
        delydot
        delzdot]
%}

tspan = [0 t]; % length of trajectory flight
Hop2.dv = Hop2.dv0_PLUS_start_burn;
Hop2.state = [Hop2.dr;Hop2.dv;Football1.rECI_target_data(end,1:3)';Football1.vECI_target_data(end,1:3)'];
[Hop2.timenew,Hop2.statenew] = ode45(@linearizedEOMs_std,tspan,Hop2.state,options,h.target,mu);

% Extract data after ODE
Hop2.rECI_target_data = [Hop2.statenew(:,7),Hop2.statenew(:,8) Hop2.statenew(:,9)];
Hop2.vECI_target_data = [Hop2.statenew(:,10),Hop2.statenew(:,11) Hop2.statenew(:,12)];
Hop2.relativePosition = [Hop2.statenew(:,1),Hop2.statenew(:,2),Hop2.statenew(:,3)]; % since z is zero whole time
Hop2.relativeVelocity = [Hop2.statenew(:,4),Hop2.statenew(:,5),Hop2.statenew(:,6)]; 

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Hop2.relativePosition(:,2),Hop2.relativePosition(:,1),'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(Hop2.relativePosition(end,2),Hop2.relativePosition(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 100 km to 40 km hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% Plot trajectory: 40km to 1km hop in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Hop2.timenew)
    hop2QXx = QXx_from_rv_ECI(Hop2.rECI_target_data(i,:)',Hop2.vECI_target_data(i,:)');
    Hop2.rECI_chaser_data(i,:) = (hop2QXx' * Hop2.relativePosition(i,:)') + Hop2.rECI_target_data(i,:)';
    Hop2.vECI_chaser_data(i,:) = (hop2QXx' * Hop2.relativeVelocity(i,:)') + Hop2.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(Hop2.rECI_target_data(:,1),Hop2.rECI_target_data(:,2),Hop2.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(Hop2.rECI_target_data(end,1),Hop2.rECI_target_data(end,2),Hop2.rECI_target_data(end,3),'*','LineWidth',5);
% p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(Hop2.rECI_chaser_data(:,1),Hop2.rECI_chaser_data(:,2),Hop2.rECI_chaser_data(:,3),'k','LineWidth',1.5,'LineStyle','--');
p4 = plot3(Hop2.rECI_chaser_data(end,1),Hop2.rECI_chaser_data(end,2),Hop2.rECI_chaser_data(end,3),'square','LineWidth',5);

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('Spacecrafts A and B at Mission t0','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' orbit','Target','Chaser', 'interpreter','latex','Location', 'best')

disp("Delta-V for Hold 1 is: " + delta_v_intoFB1 + " m/s")

missionDV.football = missionDV.hop1 + delta_v_intoFB1; % total mission DV to date; m/s
 
disp("Total Mission Delta-V so far is: " + missionDV.football + " m/s")

% Keep track of mission time
football_time = T.target / 3600; % sec to hours
MET.football = MET.hop1 + football_time; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.football + " hours")
disp("MET = " + MET.football/24 + " days")

timewarp = 86400 * (2 - MET.football/24); % SECONDS left to make even 2 days MET so far

% Let transfer trajectory be some scalar times this value. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SECOND HOP - from 40 km to 1 km

% Choose trajectory travel time.
t = timewarp; % this puts us at an even 2.0 days Mission time so far.

% Where you want to end up
Hop2.drf = [0;1;0]; % km

% relative motion for two impulse inputs
rA = Football1.rECI_target_data(end,1:3)';
rB = Football1.rECI_chaser_data(end,1:3)';
vA = Football1.vECI_target_data(end,1:3)';
vB = Football1.vECI_chaser_data(end,1:3)';

[Hop2.r_relx0, Hop2.v_relx0, Hop2.a_relx0] = rva_relative(rA,vA,rB,vB); 
% note that velocity component does not include dV for Football exit impulse

% Current dr and dv
Hop2.dr = Football1.relativePosition_Vect(end,1:3)';
Hop2.dv0 = Football1.relativeVelocity_Vect(end,1:3)'; 

% Call function to find instantaneous dv burn (start of trajectory)
%[hop1.dv0_plus,hop1.DV_0,hop1.dV_F,hop1.DV_total] = cw_twoimpulse(dr,drf,dv0,period,t);

% Call function to find instantaneous dv burn (start of trajectory)
[Hop2.dv0_PLUS_start_burn,Hop2.DV_total,Hop2.dv_F,Hop2.DV_total] = cw_twoimpulse(Hop2.dr,Hop2.drf,Hop2.dv0,period,t);

%% Plot trajectory: 40km to 1km hop in LVLH

%{
Use linearized EOMs to hop. 
state = [delx
        dely
        delz
        delxdot
        delydot
        delzdot]
%}

tspan = [0 t]; % length of trajectory flight
Hop2.dv = Hop2.dv0_PLUS_start_burn;
Hop2.state = [Hop2.dr;Hop2.dv;Football1.rECI_target_data(end,1:3)';Football1.vECI_target_data(end,1:3)'];
[Hop2.timenew,Hop2.statenew] = ode45(@linearizedEOMs_std,tspan,Hop2.state,options,h.target,mu);

% Extract data after ODE
Hop2.rECI_target_data = [Hop2.statenew(:,7),Hop2.statenew(:,8), Hop2.statenew(:,9)];
Hop2.vECI_target_data = [Hop2.statenew(:,10),Hop2.statenew(:,11), Hop2.statenew(:,12)];
Hop2.relativePosition_Vect = [Hop2.statenew(:,1),Hop2.statenew(:,2),Hop2.statenew(:,3)]; % since z is zero whole time
Hop2.relativeVelocity_Vect = [Hop2.statenew(:,4),Hop2.statenew(:,5),Hop2.statenew(:,6)]; 

relativePosition = Hop2.relativePosition_Vect(end,1:3);
Hop2.verifyDist = norm(relativePosition); % should be 1.0 km

disp("* CHECK *  Relative distance between target and chaser after maneuver is: " + Hop2.verifyDist + " km")

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Hop2.relativePosition_Vect(:,2),Hop2.relativePosition_Vect(:,1),'r','LineWidth',2)

% Chaser position after hop
p1 = plot(Hop2.relativePosition_Vect(end,2),Hop2.relativePosition_Vect(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 40 km to 1 km hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','40 km to 1 km hop', 'Chaser final position','interpreter','latex','Location', 'best')

%% Plot trajectory: 40km to 1km hop in ECI

% Convert LVLH state chaser data to ECI for plotting
% 
for i = 1:length(Hop2.timenew)
    hop2QXx = QXx_from_rv_ECI(Hop2.rECI_target_data(i,:)',Hop2.vECI_target_data(i,:)');
    Hop2.rECI_chaser_data(i,:) = (hop2QXx' * Hop2.relativePosition_Vect(i,:)') + Hop2.rECI_target_data(i,:)';
    Hop2.vECI_chaser_data(i,:) = (hop2QXx' * Hop2.relativeVelocity_Vect(i,:)') + Hop2.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET orbit and current position after 40km - 1 km maneuver
p1 = plot3(Hop2.rECI_target_data(:,1),Hop2.rECI_target_data(:,2),Hop2.rECI_target_data(:,3),'--k','LineWidth',1);
% Target position
p2 = plot3(Hop2.rECI_target_data(end,1),Hop2.rECI_target_data(end,2),Hop2.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER orbit and current position after 40km - 1 km maneuver
p3 = plot3(Hop2.rECI_chaser_data(:,1),Hop2.rECI_chaser_data(:,2),Hop2.rECI_chaser_data(:,3),'r','LineWidth',1);
p4 = plot3(Hop2.rECI_chaser_data(end,1),Hop2.rECI_chaser_data(end,2),Hop2.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('40km to 1km Hop Maneuver in ECI','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' Target orbit','Target','Chaser Orbit','Chaser', 'interpreter','latex','Location', 'best')

%% Plot all maneuvers so far in LVLH 
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
xline(0)
yline(0)

% Hop 1 trajectory
plot(relativePosition_Vect(:,2),relativePosition_Vect(:,1),'--k','LineWidth',1)

% Plot Football trajectory
plot(Football1.relativePosition_Vect(:,2),Football1.relativePosition_Vect(:,1),'--b','LineWidth',2)

% Hop trajectory
plot(Hop2.relativePosition_Vect(:,2),Hop2.relativePosition_Vect(:,1),'r','LineWidth',2)

% Chaser position after hop
p1 = plot(Hop2.relativePosition_Vect(end,2),Hop2.relativePosition_Vect(end,1),'x','LineWidth',2);
p1.Color = 'k';

% Graph pretty 
ylim([-50 50]) 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: First three maneuvers','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','','','100 km to 40 km hop','Football hold', '40 km to 1 km hop','Chaser position after hop','interpreter','latex','Location', 'best')


%<<<<<<< HEAD

% =======
%% HOP 2: 40 km to 1 km SUMMARY 
% Written by JS 10/21/23

disp(" ")
disp("--------- Hop 2: 40 km to 1 km ---------") 

disp("Delta-V for Hop 2 is: " + Hop2.DV_total + " m/s")

missionDV.hop2 = missionDV.football + Hop2.DV_total;

disp("Total Mission Delta-V so far is: " + missionDV.hop2 + " m/s")

% Keep track of mission time
currentManeuver = t/3600; % sec to hours
MET.hop2 = MET.football + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.hop2 + " hours")
disp("MET = " + MET.hop2/24 + " days")

%% HOLD 2: V-BAR station keeping at: dr = 1 km; dv = 0 m/s
% Written by JS 10/21/23

% Need to thrust continuously to stay in v-bar station keeping orbit
holdLength = 1 * 86400; % days -> seconds

% Prepare for continuous thrust; v-bar station keeping
dr = [Hop2.relativePosition_Vect(end,1);Hop2.relativePosition_Vect(end,2);0]; % dr after 40 km - 1 km hop
dv = [Hop2.relativeVelocity_Vect(end,1);Hop2.relativeVelocity_Vect(end,2);0]; % dv after 40 km - 1 km hop

%{
dv_desired = 0;
% Continuous fire for v-bar station keeping
for i = 1:t:holdLength
    [hold2.drNew,hold2.dvNew] = clohessy_wiltshire_eqs(state0,period,t);
    delta_v = dv_desired + norm(dv - hold2.dvNew); % just want the last one, thus no (i)
end

% This is the dv needed to stay at zero veloc in v-bar station keeping
hold2.dv_total = delta_v; % to report
disp("DV needed for v-bar station keeping is: " + hold2.dv_vbar/1000 + " m/s")
%}

% Another method of finding continuous firing dv
t = holdLength;
period = T.target;
drf = dr;
dv0 = dv;
[dv0_plus,DV_0,DV_f,hold2.DV_total] = cw_twoimpulse(dr,drf,dv0,period,t);

disp("DV total for the v-bar station keeping burn is: " + hold2.DV_total*1000 + " m/s"); % for hold 2

% PERFORM V-BAR STATION KEEPING BURN
tspan = [0 holdLength]; % length of trajectory flight
dv = dv0_plus;
% dr is known
rECI = Hop2.rECI_chaser_data(end,1:3);
vECI = Hop2.vECI_chaser_data(end,1:3);
hold2.state = [dr;dv;rECI';vECI'];
[~,hold2.statenew] = ode45(@linearizedEOMs_std,tspan,hold2.state,options,h.target,mu);

% Extract data after ODE
hold2.rECI_target_data = [hold2.statenew(:,7),hold2.statenew(:,8), hold2.statenew(:,9)];
hold2.vECI_target_data = [hold2.statenew(:,10),hold2.statenew(:,11), hold2.statenew(:,12)];
hold2.relativePosition_Vect = [hold2.statenew(:,1),hold2.statenew(:,2),hold2.statenew(:,3)]; % since z is zero whole time
hold2.relativeVelocity_Vect = [hold2.statenew(:,4),hold2.statenew(:,5),hold2.statenew(:,6)]; 

relativePosition = hold2.relativePosition_Vect(end,1:3);
relativeVeloc = hold2.relativeVelocity_Vect(end,1:3);
hold2.verifyDist = norm(relativePosition);      % should be 1.0 km
hold2.verifyVeloc =norm(relativeVeloc);         % should be 0 m/s

% dr should be [1 0 0]
disp("* CHECK *  Relative distance between target and chaser after maneuver is: " + hold2.verifyDist + " km")
% dv should be [0 0 0]
disp("* CHECK *  Relative velocity between target and chaser after maneuver is: " + hold2.verifyVeloc + " m/s")

% Convert LVLH state chaser data to ECI for plotting

for i = 1:length(hold2.statenew)
    hold2QXx = QXx_from_rv_ECI(hold2.rECI_target_data(i,:)',hold2.vECI_target_data(i,:)');
    hold2.rECI_chaser_data(i,:) = (hold2QXx' * hold2.relativePosition_Vect(i,:)') + hold2.rECI_target_data(i,:)';
    hold2.vECI_chaser_data(i,:) = (hold2QXx' * hold2.relativeVelocity_Vect(i,:)') + hold2.vECI_target_data(i,:)';
end

%% Plot final position of both s/c after the hold period (LVLH)
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on

% Chaser position after hold (should be 1.0 km dr; 0 dv)
p1 = plot(hold2.relativePosition_Vect(end,2),hold2.relativeVelocity_Vect(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim([-1 1]) 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: v-bar station keeping at 1 km relative distance','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Chaser relative position during hold', '','interpreter','latex','Location', 'best')

%% Plot full mission so far in LVLH

% this doesn't look any different than the last one, thankfully.
% no point in plotting.

%% Plot final position of both s/c after the hold period (ECI)

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET orbit during 1 km HOLD (hold 2)
p1 = plot3(hold2.rECI_target_data(:,1),hold2.rECI_target_data(:,2),hold2.rECI_target_data(:,3),'--r','LineWidth',1);
% Target position
p2 = plot3(hold2.rECI_target_data(end,1),hold2.rECI_target_data(end,2),hold2.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER position after 1 km v-bar station keeping hold
p3 = plot3(hold2.rECI_chaser_data(:,1),hold2.rECI_chaser_data(:,2),hold2.rECI_chaser_data(:,3),'r','LineWidth',1);
p4 = plot3(hold2.rECI_chaser_data(end,1),hold2.rECI_chaser_data(end,2),hold2.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x [km]','Interpreter','latex'); 
yLab = ylabel('y [km]','Interpreter','latex'); 
zLab = zlabel('z [km]','Interpreter','latex'); 
plotTitle = title('V-bar station keeping at 1 km in ECI','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' Target Orbit','Target','Chaser Orbit','Chaser', 'interpreter','latex','Location', 'best')

%% V-Bar Hold @ 1km: SUMMARY 
% Written by JS 10/21/23

disp(" ")
disp("--------- V-Bar station keeping at 1 km relative distance ---------") 

disp("Delta-V for this v-bar station keeping is: " + hold2.DV_total*1000 + " m/s")

missionDV.hold2 = missionDV.hop2 + hold2.DV_total*1000;

disp("Total Mission Delta-V so far is: " + missionDV.hold2 + " m/s")

% Keep track of mission time
currentManeuver = holdLength/3600; % sec to hours
MET.hold2 = MET.hop2 + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.hold2 + " hours")
disp("MET = " + MET.hold2/24 + " days")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOP 3 - From 1 km to 300 m (.30 km)
% Choose trajectory travel time.
t = 86400; % one day

% Where you want to end up
Hop3.drf = [0;0.30;0]; % km

% relative motion for two impulse inputs
rA = hold2.rECI_target_data(end,1:3)';
vA = hold2.vECI_target_data(end,1:3)';
rB = hold2.rECI_chaser_data(end,1:3)';
vB = hold2.vECI_chaser_data(end,1:3)';

[Hop3.r_relx, Hop3.v_relx, ~] = rva_relative(rA,vA,rB,vB);

% What the current relative distance and velocity (LVLH) is
Hop3.dr = Hop3.r_relx;
Hop3.dv0 = Hop3.v_relx;
% Call function to find instantaneous dv burn (start of trajectory)
[dv0_plus,DV_0,DV_f,Hop3.DV_total] = cw_twoimpulse(Hop3.dr,Hop3.drf,Hop3.dv0,period,t);

disp(" ")
disp(" --------- Hop 3: 1 km to 300m ---------") 

% Check
disp("Delta V for this hop (3): ")
disp(Hop3.DV_total/1000 + " m/s")

%% Perform maneuver: 1 km to 300 m hop in LVLH

tspan = [0 t]; % length of trajectory flight
Hop3.dv = dv0_plus;
Hop3.state = [Hop3.dr;Hop3.dv;rA;vA];
[~,Hop3.statenew] = ode45(@linearizedEOMs_std,tspan,Hop3.state,options,h.target,mu);

% Extract data after ODE
Hop3.rECI_target_data = [Hop3.statenew(:,7),Hop3.statenew(:,8) ,Hop3.statenew(:,9)];
Hop3.vECI_target_data = [Hop3.statenew(:,10),Hop3.statenew(:,11), Hop3.statenew(:,12)];
Hop3.relativePosition = [Hop3.statenew(:,1),Hop3.statenew(:,2),Hop3.statenew(:,3)]; % since z is zero whole time
Hop3.relativeVelocity = [Hop3.statenew(:,4),Hop3.statenew(:,5),Hop3.statenew(:,6)]; 

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Hop3.relativePosition(:,2)*1000,Hop3.relativeVelocity(:,1)*1000,'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(Hop3.relativePosition(end,2)*1000,Hop3.relativeVelocity(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim([-1 1])
xlim padded 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 1 km to 300 m hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% Plot trajectory: 1km to 300m hop in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Hop3.statenew)
    hop3QXx = QXx_from_rv_ECI(Hop3.rECI_target_data(i,:)',Hop3.vECI_target_data(i,:)');
    Hop3.rECI_chaser_data(i,:) = (hop3QXx' * Hop3.relativePosition(i,:)') + Hop3.rECI_target_data(i,:)';
    Hop3.vECI_chaser_data(i,:) = (hop3QXx' * Hop3.relativeVelocity(i,:)') + Hop3.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET after 1 km to 30m hop complete
p1 = plot3(Hop3.rECI_target_data(:,1),Hop3.rECI_target_data(:,2),Hop3.rECI_target_data(:,3),'--r','LineWidth',1);
p2 = plot3(Hop3.rECI_target_data(end,1),Hop3.rECI_target_data(end,2),Hop3.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER after 1 km to 30m hop complete
p3 = plot3(Hop3.rECI_chaser_data(:,1),Hop3.rECI_chaser_data(:,2),Hop3.rECI_chaser_data(:,3),'r','LineWidth',1);
p4 = plot3(Hop3.rECI_chaser_data(end,1),Hop3.rECI_chaser_data(end,2),Hop3.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('ECI frame: 1 km to 300 m hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' Target orbit','Target','Chaser orbit','Chaser', 'interpreter','latex','Location', 'best')

%% 1 km to 300 m hop: SUMMARY 
% Written by JS 10/21/23

disp(" ")
disp("--------- 1 km to 300 m hop maneuver ---------") 

disp("Delta-V for this maneuver is: " + Hop3.DV_total*1000 + " m/s")

missionDV.hop3 = missionDV.hold2 + Hop3.DV_total*1000; % next time call each phase by a seq #

disp("Total Mission Delta-V so far is: " + missionDV.hop3 + " m/s")

% Keep track of mission time
currentManeuver = t/3600; % sec to hours
MET.hop3 = MET.hold2 + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.hop3 + " hours")
disp("MET = " + MET.hop3/24 + " days")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% HOLD 3: V-BAR station keeping at: dr = 300 m; dv = 0 m/s "hold300"
% Written by JS 10/22/23

% Need to thrust continuously to stay in v-bar station keeping orbit
holdLength = 1 * 86400; % days -> seconds

% Current position and velocity
% Prepare for continuous thrust; v-bar station keeping
dr = [Hop3.relativePosition(end,1);Hop3.relativePosition(end,2);0]; % dr after 1 km to 300 m hop
dv = [Hop3.relativeVelocity(end,1);Hop3.relativeVelocity(end,2);0]; % dv after 1 km to 300 m hop
state0 = [dr;dv];


% Calculate dv necessary to keep v-bar "station keeping" orbit
t = holdLength;
period = T.target;
drf = dr;
dv0 = dv;
[dv0_plus,DV_0,DV_f,hold300.DV_total] = cw_twoimpulse(dr,drf,dv0,period,t);

disp("DV total for the v-bar station keeping burn is: " + hold300.DV_total*1000 + " m/s"); % for 300m hold

% PERFORM V-BAR STATION KEEPING BURN
tspan = [0 holdLength]; % length of trajectory flight
dv = dv0_plus;
% dr is known
rECI = Hop3.rECI_chaser_data(end,1:3); % last position of target in ECI
vECI = Hop3.vECI_chaser_data(end,1:3);
hold300.state = [dr;dv;rECI';vECI'];
[~,hold300.statenew] = ode45(@linearizedEOMs_std,tspan,hold300.state,options,h.target,mu);

% Extract data after ODE
hold300.rECI_target_data = [hold300.statenew(:,7),hold300.statenew(:,8), hold300.statenew(:,9)];
hold300.vECI_target_data = [hold300.statenew(:,10),hold300.statenew(:,11), hold300.statenew(:,12)];
hold300.relativePosition_Vect = [hold300.statenew(:,1),hold300.statenew(:,2),hold300.statenew(:,3)]; % since z is zero whole time
hold300.relativeVelocity_Vect = [hold300.statenew(:,4),hold300.statenew(:,5),hold300.statenew(:,6)]; 

relativePosition = hold300.relativePosition_Vect(end,1:3);
relativeVeloc = hold300.relativeVelocity_Vect(end,1:3);
hold300.verifyDist = norm(relativePosition);      % should be 1.0 km
hold300.verifyVeloc =norm(relativeVeloc);         % should be 0 m/s

% dr should be [1 0 0]
disp("* CHECK *  Relative distance between target and chaser after 300 m hold is: " + hold300.verifyDist + " km")
% dv should be [0 0 0]
disp("* CHECK *  Relative velocity between target and chaser after 300 m hold is: " + hold300.verifyVeloc + " m/s")

% Convert LVLH state chaser data to ECI for plotting

for i = 1:length(hold300.statenew)
    hold300QXx = QXx_from_rv_ECI(hold300.rECI_target_data(i,:)',hold300.vECI_target_data(i,:)');
    hold300.rECI_chaser_data(i,:) = (hold300QXx' * hold300.relativePosition_Vect(i,:)') + hold300.rECI_target_data(i,:)';
    hold300.vECI_chaser_data(i,:) = (hold300QXx' * hold300.relativeVelocity_Vect(i,:)') + hold300.vECI_target_data(i,:)';
end

%% Plot final position of both s/c after the hold period (LVLH)
% Written by JS 10/22/23
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on

% Chaser position after hold (should be 300 m dr; 0 dv)
p1 = plot(hold300.relativePosition_Vect(end,2)*1000,hold300.relativeVelocity_Vect(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)
yline(-150,'r')

% Error box zoomed out
%dim = [.855 .501 .03 .03]; % [x y w h]
%annotation('rectangle',dim,'Color','red')

% Error box zoomed in [275 300] m window
dim = [.75 .415 .155 .205]; % [x y w h]
annotation('rectangle',dim,'Color','red')

% Graph pretty 
ylim([-20 20]) 
xlim ([275 300]) 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: v-bar station keeping at 300 m relative distance','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Chaser relative position during hold', '','','Error box: $\pm 5$ m ','interpreter','latex','Location', 'best')
%% Plot final position of both s/c after the hold period (ECI)
% Written by JS 10/22/23

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET orbit during 300 m HOLD
p1 = plot3(hold300.rECI_target_data(:,1),hold300.rECI_target_data(:,2),hold300.rECI_target_data(:,3),'--r','LineWidth',1);
% Target position
p2 = plot3(hold300.rECI_target_data(end,1),hold300.rECI_target_data(end,2),hold300.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER position after 300 m v-bar station keeping hold
p3 = plot3(hold300.rECI_chaser_data(:,1),hold300.rECI_chaser_data(:,2),hold300.rECI_chaser_data(:,3),'r','LineWidth',1);
p4 = plot3(hold300.rECI_chaser_data(end,1),hold300.rECI_chaser_data(end,2),hold300.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('V-bar station keeping at 300 m in ECI','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' Target Orbit','Target','Chaser Orbit','Chaser', 'interpreter','latex','Location', 'best')

%% V-Bar Station Keeping Hold @ 300 m: SUMMARY 
% Written by JS 10/22/23

disp(" ")
disp("--------- V-Bar station keeping at 300 m relative distance ---------") 

disp("Delta-V for this v-bar station keeping is: " + hold300.DV_total*1000 + " m/s")

missionDV.hold300 = missionDV.hop3 + hold300.DV_total*1000;

disp("Total Mission Delta-V so far is: " + missionDV.hold300 + " m/s")

% Keep track of mission time
currentManeuver = holdLength/3600; % sec to hours
MET.hold300 = MET.hop3 + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.hold300 + " hours")
disp("MET = " + MET.hold300/24 + " days")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HOP 4 - From 300 m to 21 m
% Hop to 21 meters to prepare for the 19-21 m football hold need small dv

% Choose trajectory travel time.
t = holdLength; % one day in seconds (86400 per normal)

% Where you want to end up
Hop4.drf = [0;0.02025;0]; % km % THIS IS MODIFIED TO STAY WITHIN +- ONE METER OF 20 M for FB#2

% relative motion for two impulse inputs (last position/veloc)
rA = hold300.rECI_target_data(end,1:3)';
vA = hold300.rECI_target_data(end,1:3)';
rB = hold300.rECI_chaser_data(end,1:3)';
vB = hold300.rECI_chaser_data(end,1:3)';

[Hop4.r_relx, Hop4.v_relx, ~] = rva_relative(rA,vA,rB,vB);

% What the current relative distance and velocity (LVLH) is
Hop4.dr = hold300.relativePosition_Vect(end,1:3)';
Hop4.dv0 = hold300.relativeVelocity_Vect(end,1:3)'; 
% Call function to find instantaneous dv burn (start of trajectory)
[Hop4.dv0_plus,DV_0,DV_f,Hop4.DV_total] = cw_twoimpulse(Hop4.dr,Hop4.drf,Hop4.dv0,period,t);

disp(" ")
disp("--------- Hop 4: 300 m  to 21 m ---------") 

% Display maneuver dv:
disp("DV for the 300 m to 21 m hop is: " + Hop4.DV_total*1000 + " m/s")

%% Plot trajectory: 300 m to 21 m hop in LVLH

%{
Use linearized EOMs to hop. 
state = [delx
        dely
        delz
        delxdot
        delydot
        delzdot]
%}

tspan = [0 t]; % length of trajectory flight
Hop4.dv = Hop4.dv0_plus;
Hop4.state = [Hop4.dr;Hop4.dv;hold300.rECI_target_data(end,1:3)';hold300.vECI_target_data(end,1:3)'];
[~,Hop4.statenew] = ode45(@linearizedEOMs_std,tspan,Hop4.state,options,h.target,mu);

% Extract data after ODE
Hop4.rECI_target_data = [Hop4.statenew(:,7),Hop4.statenew(:,8), Hop4.statenew(:,9)];
Hop4.vECI_target_data = [Hop4.statenew(:,10),Hop4.statenew(:,11), Hop4.statenew(:,12)];
Hop4.relativePosition = [Hop4.statenew(:,1),Hop4.statenew(:,2),Hop4.statenew(:,3)]; % since z is zero whole time
Hop4.relativeVelocity = [Hop4.statenew(:,4),Hop4.statenew(:,5),Hop4.statenew(:,6)]; 

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Hop4.relativePosition(:,2)*1000,Hop4.relativePosition(:,1)*1000,'LineWidth',2)

% Chaser position after hop
% Plot
p1 = plot(Hop4.relativePosition(end,2)*1000,Hop4.relativePosition(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim ([-50 100]) 
xlim padded 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 300 m to 21 m hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% Plot trajectory: 300 m to 21 m hop in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Hop4.statenew)
    hop4QXx = QXx_from_rv_ECI(Hop4.rECI_target_data(i,:)',Hop4.vECI_target_data(i,:)');
    Hop4.rECI_chaser_data(i,:) = (hop4QXx' * Hop4.relativePosition(i,:)') + Hop4.rECI_target_data(i,:)';
    Hop4.vECI_chaser_data(i,:) = (hop4QXx' * Hop4.relativeVelocity(i,:)') + Hop4.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(Hop4.rECI_target_data(:,1),Hop4.rECI_target_data(:,2),Hop4.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(Hop4.rECI_target_data(end,1),Hop4.rECI_target_data(end,2),Hop4.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(Hop4.rECI_chaser_data(:,1),Hop4.rECI_chaser_data(:,2),Hop4.rECI_chaser_data(:,3),'--r','LineWidth',1);
p4 = plot3(Hop4.rECI_chaser_data(end,1),Hop4.rECI_chaser_data(end,2),Hop4.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x [km]','Interpreter','latex'); 
yLab = ylabel('y [km]','Interpreter','latex'); 
zLab = zlabel('z [km]','Interpreter','latex'); 
plotTitle = title('ECI frame: 300 m to 21 m hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' Target orbit','Target','Chaser orbit','Chaser', 'interpreter','latex','Location', 'best')

%% 300 m to 21 m hop SUMMARY 
% Written by JS 10/22/23

disp(" ")
disp("--------- 300 m to 21 m hop maneuver SUMMARY ---------") 

disp("Delta-V for this maneuver was: " + Hop4.DV_total*1000 + " m/s")

missionDV.Hop4 = missionDV.hold300 + Hop4.DV_total*1000;

disp("Total Mission Delta-V so far is: " + missionDV.Hop4 + " m/s")

% Keep track of mission time
currentManeuver = holdLength/3600; % sec to hours
MET.Hop4 = MET.hold300 + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.Hop4 + " hours")
disp("MET = " + MET.Hop4/24 + " days")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Correction burn: Need to be at true 21 m (v) and 0 (r) in LVLH

% We are not at exactly 21 m (v) and 0 (r); fix that.

% Choose trajectory travel time.
t = 3600; % one hour long correction burn

% Where you want to end up
correction.drf = [0;0.02025;0]; % km % THIS IS MODIFIED TO STAY WITHIN +- ONE METER OF 20 M for FB#2

% relative motion for two impulse inputs (last position/veloc)
rA = Hop4.rECI_target_data(end,1:3)';
vA = Hop4.rECI_target_data(end,1:3)';
rB = Hop4.rECI_chaser_data(end,1:3)';
vB = Hop4.rECI_chaser_data(end,1:3)';

[correction.r_relx, correction.v_relx, ~] = rva_relative(rA,vA,rB,vB);

% What the current relative distance and velocity (LVLH) is
correction.dr = Hop4.relativePosition(end,1:3)';
correction.dv0 = Hop4.relativeVelocity(end,1:3)'; 
% Call function to find instantaneous dv burn (start of trajectory)
[correction.dv0_plus,DV_0,DV_f,correction.DV_total] = cw_twoimpulse(correction.dr,correction.drf,correction.dv0,period,t);

disp(" ")
disp("--------- Correction Burn: to true 21 m ---------") 

% Display maneuver dv:
disp("DV for the correction burn is: " + correction.DV_total*1000 + " m/s")

% % % % % Plot trajectory: correction burn

tspan = [0 t]; % length of trajectory flight
correction.dv = correction.dv0_plus;
correction.state = [correction.dr;correction.dv;Hop4.rECI_target_data(end,1:3)';Hop4.vECI_target_data(end,1:3)'];
[~,correction.statenew] = ode45(@linearizedEOMs_std,tspan,correction.state,options,h.target,mu);

% Extract data after ODE
correction.rECI_target_data = [correction.statenew(:,7),correction.statenew(:,8), correction.statenew(:,9)];
correction.vECI_target_data = [correction.statenew(:,10),correction.statenew(:,11), correction.statenew(:,12)];
correction.relativePosition = [correction.statenew(:,1),correction.statenew(:,2),correction.statenew(:,3)]; % since z is zero whole time
correction.relativeVelocity = [correction.statenew(:,4),correction.statenew(:,5),correction.statenew(:,6)]; 

% % % Correct the ECI
for i = 1:length(correction.statenew)
    correctionQXx = QXx_from_rv_ECI(correction.rECI_target_data(i,:)',correction.vECI_target_data(i,:)');
    correction.rECI_chaser_data(i,:) = (correctionQXx' * correction.relativePosition(i,:)') + correction.rECI_target_data(i,:)';
    correction.vECI_chaser_data(i,:) = (correctionQXx' * correction.relativeVelocity(i,:)') + correction.vECI_target_data(i,:)';
end
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(correction.relativePosition(:,2)*1000,correction.relativePosition(:,1)*1000,'LineWidth',2)

% Chaser position after hop
% Plot
p1 = plot(correction.relativePosition(end,2)*1000,correction.relativePosition(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Correction hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')



%% Correction burn SUMMARY
% Written by JS 10/22/23

disp(" ")
disp("--------- Correction burn at ~21 meters ---------") 

disp("Delta-V for this maneuver was: " + correction.DV_total*1000 + " m/s")

missionDV.correction = missionDV.Hop4 + correction.DV_total*1000;

disp("Total Mission Delta-V so far is: " + missionDV.correction + " m/s")

% Keep track of mission time
currentManeuver = t/3600; % sec to hours
MET.correction = MET.Hop4 + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.correction + " hours")
disp("MET = " + MET.correction/24 + " days")
disp("We will need an additional 23 hours to make an even 7 days. ")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time is off after the correction burn. Let's get back on track with whole
% numbers (we will take up the slack with the final burn)

% Need additional 23 hours. ***
%% Football between 19 and 21

% Set major axis of football orbit (between 19-21 meters away with dv)
b = 0.1/1000; % one meter football hold around 20 m 
Football2.xdot0  = ((b)) * n.Target; % altitude direction
dv_FootBall2_LVLH = [Football2.xdot0 0 0];
Football2.dvChaser0 = correction.relativeVelocity(end,1:3) + dv_FootBall2_LVLH;

disp(" ")
disp("------------------")
disp("Delta V from Football Manuever 2: " + norm(Football2.dvChaser0)*1000 + " m/s")

Football2.dr_t = correction.relativePosition(end,1:3);
Football2.dv_t = dv_FootBall2_LVLH;
Football2.dr0 = Football2.dr_t;
Football2.dv0 = Football2.dv_t;

t = T.target; % 23 hours (to make up for that 1 hour correction burn earlier)

Football2.relativePosition = zeros(1,3);
Football2.relativeVelocity = zeros(1,3);

% Inputs for relative position calculation
rA = correction.rECI_target_data(end,1:3)';
vA = correction.vECI_target_data(end,1:3)';
rB = correction.rECI_chaser_data(end,1:3)';
vB = correction.vECI_chaser_data(end,1:3)';

% relative motion
[Football2.r_relx, Football2.v_relx, ~] = rva_relative(rA,vA,rB,vB);

tspan = [0 t]; % length of trajectory flight
Football2.dr = Football2.relativePosition(end,1:3);
Football2.dv = Football2.relativePosition(end,1:3);

state_FB2 = [Football2.dr0';Football2.dv0';...
correction.rECI_target_data(end,1:3)'; correction.vECI_target_data(end,1:3)'];
[~,Football2.statenew] = ode45(@linearizedEOMs_std,tspan,state_FB2,options,h.target,mu);

% Extract data after ODE
Football2.rECI_target_data = [Football2.statenew(:,7),Football2.statenew(:,8), Football2.statenew(:,9)];
Football2.vECI_target_data = [Football2.statenew(:,10),Football2.statenew(:,11), Football2.statenew(:,12)];
Football2.relativePosition = [Football2.statenew(:,1),Football2.statenew(:,2),Football2.statenew(:,3)]; % since z is zero whole time
Football2.relativeVelocity = [Football2.statenew(:,4),Football2.statenew(:,5),Football2.statenew(:,6)]; 

% Plot in LVLH
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Football2.relativePosition(:,2)*1000,Football2.relativePosition(:,1)*1000,'LineWidth',2)

% Chaser position after hop
% Plot
p1 = plot(Football2.relativePosition(end,2)*1000,Football2.relativePosition(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded
xlim padded 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Football between 19 m and 21 m','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')


%% NEXT: Plot football 20 m hold in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Football2.statenew)
    FB2QXx = QXx_from_rv_ECI(Football2.rECI_target_data(i,:)',Football2.vECI_target_data(i,:)');
    Football2.rECI_chaser_data(i,:) = (FB2QXx' * Football2.relativePosition(i,:)') + Football2.rECI_target_data(i,:)';
    Football2.vECI_chaser_data(i,:) = (FB2QXx' * Football2.relativeVelocity(i,:)') + Football2.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET after football 20 m hold
p1 = plot3(Football2.rECI_target_data(:,1),Football2.rECI_target_data(:,2),Football2.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(Football2.rECI_target_data(end,1),Football2.rECI_target_data(end,2),Football2.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER after football 20 m hold
p3 = plot3(Football2.rECI_chaser_data(:,1),Football2.rECI_chaser_data(:,2),Football2.rECI_chaser_data(:,3),'--r','LineWidth',1);
p4 = plot3(Football2.rECI_chaser_data(end,1),Football2.rECI_chaser_data(end,2),Football2.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('ECI frame: Football between 19 m and 21 m','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('','Target Orbit','Target','Chaser Path', 'Chaser', 'interpreter','latex','Location', 'best')

%% Football 20m hold SUMMARY

% Written by JS 10/22/23

disp(" ")
disp("--------- Football trajectory hold at ~20 m ---------") 

Football2.DV_total = norm(Football2.dvChaser0)*1000;
disp("Delta-V for this maneuver was: " + Football2.DV_total + " m/s")

missionDV.Football2 = missionDV.correction + Football2.DV_total*1000;

disp("Total Mission Delta-V so far is: " + missionDV.Football2 + " m/s")

% Keep track of mission time
currentManeuver = t/3600; % sec to hours
MET.Football2 = MET.correction + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.Football2 + " hours")
disp("MET = " + MET.Football2/24 + " days")

t2tenDays.days = 10 - MET.Football2/24; % days
t2tenDays.seconds = t2tenDays.days*24*3600; % sec
disp("We will need an additional " + t2tenDays.days + " days to make our ten day goal")
disp("Which is: " + t2tenDays.seconds + " seconds to make our ten day goal")

% Check time
equalsTen = MET.Football2/(24) + t2tenDays.days; % yes


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % NEXT MANEUVER % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Terminal maneuver: v-bar approach to rendezvous; "term"
close all; clc;
% Burn off old trajectory into super slow v-bar approach.
t = t2tenDays.seconds; % remaining time in mission (brings us to 10.000 days)
% %%%%%%%%%%%%%%%%% Find burn dv
% Where you want to end up
term.drf = [0;0;0]; % km

% Current dr, dv
term.dr = Football2.relativePosition(end,1:3)';
term.dv0 = Football2.relativeVelocity(end,1:3)'; 
term.dr(1) = 0;
% Call function to find instantaneous dv burn (start of trajectory)
[term.dv0_PLUS_start_burn,term.DV_total,term.dv_F,term.DV_total] = cw_twoimpulse(term.dr,term.drf,term.dv0,period,t);

% %%%%%%%%%%%%%%% Propogate trajectory
tspan = [0 t]; % length of trajectory flight

% Initial ECI pos, veloc
term.rECI_target = Football2.rECI_target_data(end,1:3)';
term.vECI_target = Football2.vECI_target_data(end,1:3)';

% Initial Relative position and veloc; CHASER
vc = term.dr(2)/ t; % km/s
term.dv = Football2.relativeVelocity(end,1:3)'; % km/s
term.dv(2) = -vc;
term.dv(1) = 0;

% Call v-bar propogate function
state = [term.dr;term.dv;term.rECI_target;term.vECI_target];
[~,term.statenew] = ode45(@vbar_approach2,tspan,state,options,n.Target,mu);

% Extract data after ODE

term.relativePosition = [term.statenew(:,1),term.statenew(:,2), term.statenew(:,3)];
term.relativeVelocity = [term.statenew(:,4),term.statenew(:,5), term.statenew(:,6)]; 
term.rECI_target_data = [term.statenew(:,7),term.statenew(:,8), term.statenew(:,9)];
term.vECI_target_data = [term.statenew(:,10),term.statenew(:,11), term.statenew(:,12)];

disp("Final position rel is: ")
disp(term.relativePosition(end,1:3))

disp("Final velocity rel is: ")
disp(term.relativeVelocity(end,1:3))

% % % % Plot final approach in LVLH

% Plot in LVLH
figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(term.relativePosition(:,2)*1000,term.relativePosition(:,1)*1000,'LineWidth',2)

% Chaser position after hop
% Plot
p1 = plot(term.relativePosition(end,2)*1000,term.relativePosition(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded
xlim padded 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Final approach: $v$-bar','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','$v$-bar approach', 'Rendezvous!','interpreter','latex','Location', 'best')


% % % Plot final approach in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(term.statenew)
    termQXx = QXx_from_rv_ECI(term.rECI_target_data(i,:)',term.vECI_target_data(i,:)');
    term.rECI_chaser_data(i,:) = (termQXx' * term.relativePosition(i,:)') + term.rECI_target_data(i,:)';
    term.vECI_chaser_data(i,:) = (termQXx' * term.relativeVelocity(i,:)') + term.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET after rendezvous
p1 = plot3(term.rECI_target_data(:,1),term.rECI_target_data(:,2),term.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(term.rECI_target_data(end,1),term.rECI_target_data(end,2),term.rECI_target_data(end,3),'square','LineWidth',2);
p2.Color = 'b';

% Show CHASER after rendezvous
p3 = plot3(term.rECI_chaser_data(:,1),term.rECI_chaser_data(:,2),term.rECI_chaser_data(:,3),'--r','LineWidth',1);
p4 = plot3(term.rECI_chaser_data(end,1),term.rECI_chaser_data(end,2),term.rECI_chaser_data(end,3),'x','LineWidth',2);
p4.Color = 'k';

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x [km]','Interpreter','latex'); 
yLab = ylabel('y [km]','Interpreter','latex'); 
zLab = zlabel('z [km]','Interpreter','latex'); 
plotTitle = title('ECI frame: Final $v$-bar approach and rendezvous','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('','Target Orbit','Target','Chaser Path', 'Chaser', 'interpreter','latex','Location', 'best')


disp("Delta V from Terminal Vbar Approach: " + term.DV_total*1000 + " m/s")
%% Plot last few hops of Mission in LVLH

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
xline(0)
yline(0)

% Correction Burn
plot(correction.relativePosition(:,2)*1000,correction.relativePosition(:,1)*1000,'--k','LineWidth',1)

% 20 m football hold
plot(Football2.relativePosition(:,2)*1000,Football2.relativePosition(:,1)*1000,'r','LineWidth',1)

% v-bar final approach
plot(term.relativePosition(:,2)*1000,term.relativePosition(:,1)*1000,'b','LineWidth',2)

% Chaser position at rendezvous!
p1 = plot(Football2.relativePosition(end,2)*1000,Football2.relativePosition(end,1)*1000,'x','LineWidth',2);
p1.Color = 'k';

% Graph pretty 
ylim padded
xlim padded 
xLab = xlabel('Downrange [m]','Interpreter','latex'); 
yLab = ylabel('Altitude [m]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Close approach','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','','','Position correction burn','Football hold at 20 m', 'Final approach','Chaser position after hop','interpreter','latex','Location', 'best')
 
%% V-bar approach (terminal) SUMMARY
% Written by JS 10/22/23
disp(" ")
disp("--------- Final approach (v-bar) ---------") 

disp("Delta-V for terminal v-bar maneuver was: " + term.DV_total*1000 + " m/s")

missionDV.term = missionDV.Football2 + term.DV_total*1000;

disp("Total Mission Delta-V so far is: " + missionDV.term + " m/s")

% Keep track of mission time
currentManeuver = t/3600; % sec to hours
MET.term = MET.Football2 + currentManeuver; % MISSION ELAPSED TIME (hours)
disp(" ")
disp("**** Mission Elapsed Time ****")
disp("MET = " + MET.term + " hours")
disp("MET = " + MET.term/24 + " days")


% Report final mission delta-v used (should be in the 3ish m/s range)
% >>>>>>> main
