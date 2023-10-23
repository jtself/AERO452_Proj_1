function dydt = vbar_approach2(time,state,n,mu)
%{
For use in ode45 for propogating orbits through when using a terminal or
v-bar approach for s/c rendezvous. 

ASSUMES CIRCULAR ORBITS!!

INPUTS: 
    state of initial conditions:
        state(x:y) = something
        state(and so on)...
    n = mean motion = 2pi/T
    mu = grav param

OUTPUT: dydt for ode45 to numerically integrate
%}

y = state; 

rvect = state(7:9);
v = state(10:12);
h = norm(cross(rvect,v));
R = norm(rvect);

% DEPUTY, chaser
% velocity
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);

% Linearized EOMs for acceleration IN V-BAR APPROACH SETTING 
dydt(4) = 3*n^2*y(1) + 2*n*y(5) - 2*n*y(5);
dydt(5) = 0; %-2*n*y(4) + 0;
dydt(6) = 0; %-n^2*y(3) + 0;


% CHIEF, target
% Just propogating forward like normal using acceleration function (pig function).
dydt(7) = y(10);
dydt(8) = y(11);
dydt(9) = y(12);
dydt(10) = -mu * y(7) / R^3;
dydt(11) = -mu * y(8) / R^3;
dydt(12) = -mu * y(9) / R^3;

dydt = dydt';

end % function