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
r.target = alt.target + r_earth; % km
ecc.target = 0.0000178;
theta.target = 162.2886; % deg
T.target = 2*pi*sqrt(r.target^3/mu);

h.target = findh(r.target,mu,ecc.target,theta.target);

[rECI.target,vECI.target] = r_and_v_from_COEs(RAAN.target,inc.target,w.target,h.target,ecc.target,theta.target);
h_target_vector = cross(rECI.target,vECI.target);
%% CHASER initial conditions

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
disp("Relative distance between target and chaser at mission t0 is: " + rho_rel + " km")

%% Plot mission time T0 orbit(s)
tspan = [0 T.target]; % seconds = 1 day
options = odeset('RelTol', 1e-8, 'AbsTol',1e-8);
state0.target = [rECI.target,vECI.target];

% Propogate 1 period TARGET
[newtime0.target, newstate0.target] = ode45(@coast_ODE,tspan,state0.target,options,mu);

figure(1)
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(newstate0.target(:,1),newstate0.target(:,2),newstate0.target(:,3),'r','LineWidth',2);
p2 = plot3(newstate0.target(end,1),newstate0.target(end,2),newstate0.target(end,3),'*','LineWidth',5);
p2.Color = 'k';

% Show CHASER at mission time t0
plot3(rECI.chaser(1),rECI.chaser(2),rECI.chaser(3),'*','LineWidth',5)

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
legend('','Common orbit','Target','Chaser', 'interpreter','latex','Location', 'best')

%% Plot mission t0 configuration in LVLH

vbar0.chaser = rLVLH.target(2) - rLVLH.chaser(2);
rbar0.chaser = rLVLH.target(1) - rLVLH.chaser(1);

figure(2)
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

%% FIRST RENDEZVOUS MANEUVER: 100km to 40km hop
period = T.target; 

% Choose trajectory travel time.
hours = 6; % your choice
t = 3600*hours;

% Where you want to end up
drf = [0;40;0];

% What the current relative distance and velocity (LVLH) is
dr = r_relx;
dv0 = v_relx;

% Call function to find instantaneous dv burn (start of trajectory)
[hop1.dv0_PLUS_start_burn,hop1.dvf_MINUS_off_burn,hop1.deltaV,hop1.deltaV_afterBurn] = cw_twoimpulse(dr,drf,dv0,period,t);

disp("--------- Hop: 100km to 40km ---------") 

% dv0_PLUS_start_burn % just display it
disp("Delta V: ")
disp(hop1.deltaV)
disp("Delta V After Burn: ")
disp(hop1.deltaV_afterBurn)

%% Plot trajectory: 100km to 40km hop

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
dv = hop1.dv0_PLUS_start_burn;
state = [dr;dv;rECI.target;vECI.target];
[timenew,statenew] = ode45(@linearizedEOMs_std,tspan,state,options,h.target,mu);

% Extract data after ODE
hop1.rECI_target_data = [statenew(:,7),statenew(:,8) statenew(:,9)];
hop1.vECI_target_data = [statenew(:,10),statenew(:,11) statenew(:,12)];
relativePosition = [statenew(:,1),statenew(:,2),statenew(:,3)]; % since z is zero whole time
relativeVelocity = [statenew(:,4),statenew(:,5),statenew(:,6)]; 

figure(3)
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(relativePosition(:,2),relativePosition(:,1),'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(relativePosition(end,2),relativePosition(end,1),'x','LineWidth',2);
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

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(timenew)
    hop1QXx = QXx_from_rv_ECI(hop1.rECI_target_data(i,:)',hop1.vECI_target_data(i,:)');
    hop1.rECI_chaser_data(i,:) = (hop1QXx' * relativePosition(i,:)') + hop1.rECI_target_data(i,:)';
    hop1.vECI_chaser_data(i,:) = (hop1QXx' * relativeVelocity(i,:)') + hop1.vECI_target_data(i,:)';
end

figure(4)
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(hop1.rECI_target_data(:,1),hop1.rECI_target_data(:,2),hop1.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(hop1.rECI_target_data(end,1),hop1.rECI_target_data(end,2),hop1.rECI_target_data(end,3),'*','LineWidth',5);
% p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(hop1.rECI_chaser_data(:,1),hop1.rECI_chaser_data(:,2),hop1.rECI_chaser_data(:,3),'k','LineWidth',1.5,'LineStyle','--');
p4 = plot3(hop1.rECI_chaser_data(end,1),hop1.rECI_chaser_data(end,2),hop1.rECI_chaser_data(end,3),'*','LineWidth',5);

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
legend('',' orbit','Target','Chaser', 'interpreter','latex','Location', 'best')


%% NEXT: Burn off this hop trajectory onto v-bar in a static fashion FOOTBALL

n.Target = (2*pi)/T.target; % 1/seconds
b = relativePosition(end,2)/2;
Football1.xdot0  = ((b)) * n.Target; % altitude direction
dv_FootBall1_LVLH = [Football1.xdot0 0 0];
Football1.dvChaser0 = relativeVelocity(end,1:3) + dv_FootBall1_LVLH;

disp(" ")
disp("------------------ ")
disp("Delta V from Football Manuever 1: " + norm(Football1.dvChaser0) + " km/s")

dr_t = relativePosition(end,1:3);
dv_t = dv_FootBall1_LVLH;
dr0 = dr_t;
dv0 = dv_t;

t = T.target/50000;

Football1.relativePosition = zeros(1:3);
Football1.relativeVelocity = zeros(1:3);

% ----- ALTERNATE METHOD USING CW MATRIX ---------
% for i = 1:50000
% [dr_t,dv_t] = CW_Matrix(t,dr0',dv0',n.Target);
% Football1.relativePosition(i,1:3) = dr_t;
% Football1.relativeVelocity(i,1:3) = dv_t;
% dr0 = dr_t';
% dv0 = dv_t';
% end


% relative motion
[Football1.r_relx0, Football1.v_relx0, Football1.a_relx0] = rva_relative(...
    hop1.rECI_chaser_data(end,1:3)',hop1.vECI_chaser_data(end,1:3)',...
    hop1.rECI_target_data(end,1:3)',hop1.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse


tspan = [0 T.target]; % length of trajectory flight
Football1.dr0 = relativePosition(end,1:3);
Football1.dv0 = relativeVelocity(end,1:3) + dv_FootBall1_LVLH;
dr = relativePosition(end,1:3);
dv = Football1.dvChaser0;
state_FB1 = [dr0';dv0';...
    hop1.rECI_target_data(end,1:3)';hop1.vECI_target_data(end,1:3)'];
[Football1.timenew,Football1.statenew] = ode45(@linearizedEOMs_std,tspan,state_FB1,options,h.target,mu);


% Extract data after ODE
Football1.rECI_target_data = [Football1.statenew(:,7),Football1.statenew(:,8) Football1.statenew(:,9)];
Football1.vECI_target_data = [Football1.statenew(:,10),Football1.statenew(:,11) Football1.statenew(:,12)];
Football1.relativePosition = [Football1.statenew(:,1),Football1.statenew(:,2),Football1.statenew(:,3)]; % since z is zero whole time
Football1.relativeVelocity = [Football1.statenew(:,4),Football1.statenew(:,5),Football1.statenew(:,6)]; 

figure(5)
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Football1.relativePosition(:,2),Football1.relativePosition(:,1),'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(Football1.relativePosition(end,2),Football1.relativePosition(end,1),'x','LineWidth',2);
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
legend('Target','Chaser Trajectory (football orbit)', 'Chaser final position','interpreter','latex','Location', 'best')

%% NEXT: Plot the first FOOTBALL maneuver in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Football1.timenew)
    FB1QXx = QXx_from_rv_ECI(Football1.rECI_target_data(i,:)',Football1.vECI_target_data(i,:)');
    Football1.rECI_chaser_data(i,:) = (FB1QXx' * Football1.relativePosition(i,:)') + Football1.rECI_target_data(i,:)';
    Football1.vECI_chaser_data(i,:) = (FB1QXx' * Football1.relativeVelocity(i,:)') + Football1.vECI_target_data(i,:)';
end

figure(6)
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(Football1.rECI_target_data(:,1),Football1.rECI_target_data(:,2),Football1.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(Football1.rECI_target_data(end,1),Football1.rECI_target_data(end,2),Football1.rECI_target_data(end,3),'*','LineWidth',5);
% p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(Football1.rECI_chaser_data(:,1),Football1.rECI_chaser_data(:,2),Football1.rECI_chaser_data(:,3),'k','LineWidth',1.5,'LineStyle','--');
p4 = plot3(Football1.rECI_chaser_data(end,1),Football1.rECI_chaser_data(end,2),Football1.rECI_chaser_data(end,3),'*','LineWidth',5);

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

%% SECOND HOP - from 40 km to 1 km

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
disp("--------- Hop 2: 40km to 1km ---------") 

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
plotTitle = title('LVLH frame: 40 km to 1 km hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Second Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

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
plotTitle = title('HOLD 3: 40km to 1km Hop Maneuver in ECI','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' orbit','Target','Chaser', 'interpreter','latex','Location', 'best')


%% HOP 3 - From 1 km to 30 m (.03 km)

% Choose trajectory travel time.
hours = 6; % your choice
t = 3600*hours;

% Where you want to end up
Hop3.drf = [0;30/1000;0];

% relative motion for two impulse inputs
[Hop3.r_relx0, Hop3.v_relx0, Hop3.a_relx0] = rva_relative(...
    Hop2.rECI_chaser_data(end,1:3)',Hop2.vECI_chaser_data(end,1:3)',...
    Hop2.rECI_target_data(end,1:3)',Hop2.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse


% What the current relative distance and velocity (LVLH) is
Hop3.dr = Hop2.relativePosition(end,1:3)';
Hop3.dv0 = Hop2.relativeVelocity(end,1:3)'; 
% Call function to find instantaneous dv burn (start of trajectory)
[Hop3.dv0_PLUS_start_burn,Hop3.dvf_MINUS_off_burn,Hop3.deltaV,Hop3.deltaV_afterBurn] = cw_twoimpulse(Hop3.dr,Hop3.drf,Hop3.dv0,period,t);

disp(" ")
disp(" --------- Hop 3: 30km to 0.03km (30m) ---------") 

% dv0_PLUS_start_burn % just display it
disp("Delta V: ")
disp(Hop3.deltaV)
disp("Delta V After Burn: ")
disp(Hop3.deltaV_afterBurn)

% SECOND HOP - from 40 km to 1 km

% Choose trajectory travel time.
hours = 6; % your choice
t = 3600*hours;

%% Plot trajectory: 1km to 30 meters hop in LVLH

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
Hop3.dv = Hop3.dv0_PLUS_start_burn;
Hop3.state = [Hop3.dr;Hop3.dv;Hop2.rECI_target_data(end,1:3)';Hop2.vECI_target_data(end,1:3)'];
[Hop3.timenew,Hop3.statenew] = ode45(@linearizedEOMs_std,tspan,Hop3.state,options,h.target,mu);

% Extract data after ODE
Hop3.rECI_target_data = [Hop3.statenew(:,7),Hop3.statenew(:,8) Hop3.statenew(:,9)];
Hop3.vECI_target_data = [Hop3.statenew(:,10),Hop3.statenew(:,11) Hop3.statenew(:,12)];
Hop3.relativePosition = [Hop3.statenew(:,1),Hop3.statenew(:,2),Hop3.statenew(:,3)]; % since z is zero whole time
Hop3.relativeVelocity = [Hop3.statenew(:,4),Hop3.statenew(:,5),Hop3.statenew(:,6)]; 

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Hop3.relativePosition(:,2),Hop3.relativePosition(:,1),'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(Hop3.relativePosition(end,2),Hop3.relativePosition(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 1 km to 30 m hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% Plot trajectory: 1km to 30m hop in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Hop3.timenew)
    hop3QXx = QXx_from_rv_ECI(Hop3.rECI_target_data(i,:)',Hop3.vECI_target_data(i,:)');
    Hop3.rECI_chaser_data(i,:) = (hop3QXx' * Hop3.relativePosition(i,:)') + Hop3.rECI_target_data(i,:)';
    Hop3.vECI_chaser_data(i,:) = (hop3QXx' * Hop3.relativeVelocity(i,:)') + Hop3.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(Hop3.rECI_target_data(:,1),Hop3.rECI_target_data(:,2),Hop3.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(Hop3.rECI_target_data(end,1),Hop3.rECI_target_data(end,2),Hop3.rECI_target_data(end,3),'*','LineWidth',5);
% p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(Hop3.rECI_chaser_data(:,1),Hop3.rECI_chaser_data(:,2),Hop3.rECI_chaser_data(:,3),'k','LineWidth',1.5,'LineStyle','--');
p4 = plot3(Hop3.rECI_chaser_data(end,1),Hop3.rECI_chaser_data(end,2),Hop3.rECI_chaser_data(end,3),'square','LineWidth',5);

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('ECI frame: 1 km to 30 m hop','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' orbit','Target','Chaser', 'interpreter','latex','Location', 'best')


%% HOP 4 - From 30 m to 21 m (.03 km)
% Hop to 21 Meters so we can initiate the Football Orbit After

% Choose trajectory travel time.
hours = 6; % your choice
t = 3600*hours;

% Where you want to end up
Hop4.drf = [0;21/1000;0];

% relative motion for two impulse inputs
[Hop4.r_relx0, Hop4.v_relx0, Hop4.a_relx0] = rva_relative(...
    Hop3.rECI_chaser_data(end,1:3)',Hop3.vECI_chaser_data(end,1:3)',...
    Hop3.rECI_target_data(end,1:3)',Hop3.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse


% What the current relative distance and velocity (LVLH) is
Hop4.dr = Hop3.relativePosition(end,1:3)';
Hop4.dv0 = Hop3.relativeVelocity(end,1:3)'; 
% Call function to find instantaneous dv burn (start of trajectory)
[Hop4.dv0_PLUS_start_burn,Hop4.dvf_MINUS_off_burn,Hop4.deltaV,Hop4.deltaV_afterBurn] = cw_twoimpulse(Hop4.dr,Hop4.drf,Hop4.dv0,period,t);

disp(" ")
disp("--------- Hop 4: 30km to 0.03km (30m) ---------") 

% dv0_PLUS_start_burn % just display it
disp("Delta V: ")
disp(Hop4.deltaV)
disp("Delta V After Burn: ")
disp(Hop4.deltaV_afterBurn)

% SECOND HOP - from 40 km to 1 km

% Choose trajectory travel time.
hours = 6; % your choice
t = 3600*hours;

%% Plot trajectory: 1km to 30 meters hop in LVLH

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
Hop4.dv = Hop4.dv0_PLUS_start_burn;
Hop4.state = [Hop4.dr;Hop4.dv;Hop3.rECI_target_data(end,1:3)';Hop3.vECI_target_data(end,1:3)'];
[Hop4.timenew,Hop4.statenew] = ode45(@linearizedEOMs_std,tspan,Hop4.state,options,h.target,mu);

% Extract data after ODE
Hop4.rECI_target_data = [Hop4.statenew(:,7),Hop4.statenew(:,8) Hop4.statenew(:,9)];
Hop4.vECI_target_data = [Hop4.statenew(:,10),Hop4.statenew(:,11) Hop4.statenew(:,12)];
Hop4.relativePosition = [Hop4.statenew(:,1),Hop4.statenew(:,2),Hop4.statenew(:,3)]; % since z is zero whole time
Hop4.relativeVelocity = [Hop4.statenew(:,4),Hop4.statenew(:,5),Hop4.statenew(:,6)]; 

figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Hop4.relativePosition(:,2),Hop4.relativePosition(:,1),'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(Hop4.relativePosition(end,2),Hop4.relativePosition(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: 30 m to 21 m hop','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% Plot trajectory: 1km to 30m hop in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Hop4.timenew)
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
p2 = plot3(Hop4.rECI_target_data(end,1),Hop4.rECI_target_data(end,2),Hop4.rECI_target_data(end,3),'*','LineWidth',5);
% p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(Hop4.rECI_chaser_data(:,1),Hop4.rECI_chaser_data(:,2),Hop4.rECI_chaser_data(:,3),'k','LineWidth',1.5,'LineStyle','--');
p4 = plot3(Hop4.rECI_chaser_data(end,1),Hop4.rECI_chaser_data(end,2),Hop4.rECI_chaser_data(end,3),'square','LineWidth',5);

% Graph pretty 
ylim padded 
xlim padded 
zlim padded
xLab = xlabel('x','Interpreter','latex'); 
yLab = ylabel('y','Interpreter','latex'); 
zLab = zlabel('z','Interpreter','latex'); 
plotTitle = title('ECI frame: 30 to 21 m hop','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab, zLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab, zLab],'FontSize', 14) 
grid on 
legend('',' orbit','Target','Chaser', 'interpreter','latex','Location', 'best')




%% Football between 19 and 21


b = 1/1000;
Football2.xdot0  = ((b)) * n.Target; % altitude direction
dv_FootBall2_LVLH = [Football2.xdot0 0 0];
Football2.dvChaser0 = Hop4.relativeVelocity(end,1:3) + dv_FootBall2_LVLH;

disp(" ")
disp("------------------")
disp("Delta V from Football Manuever 2: " + norm(Football2.dvChaser0) + " km/s")

Football2.dr_t = Hop4.relativePosition(end,1:3);
Football2.dv_t = dv_FootBall2_LVLH;
Football2.dr0 = Football2.dr_t;
Football2.dv0 = Football2.dv_t;

t = T.target/50000;

Football2.relativePosition = zeros(1:3);
Football2.relativeVelocity = zeros(1:3);

% % ----- ALTERNATE METHOD USING CW MATRIX ---------
%  for i = 1:50000
%  [Football2.dr_t,Football2.dv_t] = CW_Matrix(t,Football2.dr0',Football2.dv0',n.Target);
%  Football2.relativePosition(i,1:3) = Football2.dr_t;
%  Football2.relativeVelocity(i,1:3) = Football2.dv_t;
%  Football2.dr0 = Football2.dr_t';
%  Football2.dv0 = Football2.dv_t';
%  end

% 
% relative motion
[Football2.r_relx0, Football2.v_relx0, Football2.a_relx0] = rva_relative(...
    Hop4.rECI_chaser_data(end,1:3)',Hop4.vECI_chaser_data(end,1:3)',...
    Hop4.rECI_target_data(end,1:3)',Hop4.vECI_target_data(end,1:3)'); % note that velocity component does not include dV for Football entry impulse


tspan = [0 T.target]; % length of trajectory flight
Football2.dr = relativePosition(end,1:3);
Football2.dv = Football2.dvChaser0;
state_FB2 = [Football2.dr0';Football2.dv0';...
    Hop4.rECI_target_data(end,1:3)';Hop4.vECI_target_data(end,1:3)'];
[Football2.timenew,Football2.statenew] = ode45(@linearizedEOMs_std,tspan,state_FB2,options,h.target,mu);


% Extract data after ODE
Football2.rECI_target_data = [Football2.statenew(:,7),Football2.statenew(:,8) Football2.statenew(:,9)];
Football2.vECI_target_data = [Football2.statenew(:,10),Football2.statenew(:,11) Football2.statenew(:,12)];
Football2.relativePosition = [Football2.statenew(:,1),Football2.statenew(:,2),Football2.statenew(:,3)]; % since z is zero whole time
Football2.relativeVelocity = [Football2.statenew(:,4),Football2.statenew(:,5),Football2.statenew(:,6)]; 


figure()
% target, center of LVLH frame
plot(0,0,'square','Linewidth',2)
hold on
% Hop trajectory
plot(Football2.relativePosition(:,2),Football2.relativePosition(:,1),'LineWidth',2)

% Chaser position after hop

% Plot
p1 = plot(Football2.relativePosition(end,2),Football2.relativePosition(end,1),'x','LineWidth',2);
p1.Color = 'k';
xline(0)
yline(0)

% Graph pretty 
ylim padded 
xlim padded 
xLab = xlabel('Downrange [km]','Interpreter','latex'); 
yLab = ylabel('Altitude [km]','Interpreter','latex'); 
plotTitle = title('LVLH frame: Football between 19 m and 21 m','interpreter','latex'); 
set(plotTitle,'FontSize',14,'FontWeight','bold') 
set(gca,'FontName','Palatino Linotype') 
set([xLab, yLab],'FontName','Palatino Linotype') 
set(gca,'FontSize', 9) 
set([xLab, yLab],'FontSize', 14) 
grid on 
legend('Target','Hop maneuver', 'Chaser final position','interpreter','latex','Location', 'best')

%% NEXT: Plot the first FOOTBALL maneuver in ECI

% Convert LVLH state data of the chaser on the first hop to ECI
% 
for i = 1:length(Football2.timenew)
    FB2QXx = QXx_from_rv_ECI(Football2.rECI_target_data(i,:)',Football2.vECI_target_data(i,:)');
    Football2.rECI_chaser_data(i,:) = (FB2QXx' * Football2.relativePosition(i,:)') + Football2.rECI_target_data(i,:)';
    Football2.vECI_chaser_data(i,:) = (FB2QXx' * Football2.relativeVelocity(i,:)') + Football2.vECI_target_data(i,:)';
end

figure
   h1 = gca;
   earth_sphere(h1)
   hold on

% TARGET at mission start time, t0
p1 = plot3(Football2.rECI_target_data(:,1),Football2.rECI_target_data(:,2),Football2.rECI_target_data(:,3),'r','LineWidth',2);
p2 = plot3(Football2.rECI_target_data(end,1),Football2.rECI_target_data(end,2),Football2.rECI_target_data(end,3),'*','LineWidth',5);
% p2.Color = 'b';

% Show CHASER at mission time t0
p3 = plot3(Football2.rECI_chaser_data(:,1),Football2.rECI_chaser_data(:,2),Football2.rECI_chaser_data(:,3),'k','LineWidth',1.5,'LineStyle','--');
p4 = plot3(Football2.rECI_chaser_data(end,1),Football2.rECI_chaser_data(end,2),Football2.rECI_chaser_data(end,3),'*','LineWidth',5);

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




