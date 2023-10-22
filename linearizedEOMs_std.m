function [dydt] = linearizedEOMs_std(time,state,h,mu)

% This is for rendezvous calculations for "sufficiently close" approach
% using Linearized EOMS (no CW); since there is some elliptical-ness to the
% target orbit. 

% for use in ode45

%{
Note that y = state = [ [r] (3x1)
                        [v] (3x1)
                        delx
                        dely
                        delz
                        delxdot
                        delydot
                        delzdot]
SEE CURTIS PAGE 363, Equation 7.36
%}

y = state; % rename for simplicity

rvect = state(7:9);
V = state(10:12);
R = norm(rvect);

% DEPUTY, chaser
% velocity
dydt(1) = y(4);
dydt(2) = y(5);
dydt(3) = y(6);

% Linearized EOMs
dydt(4) = ((2*mu)/(R^3) + (h^2 / R^4)) * y(1) - (2 * dot(V,rvect)*(h/R^4)) * y(2) + 2*(h/R^2)*y(5);
dydt(5) = ((h^2)/(R^4) - (mu / R^3)) * y(2) + (2 * dot(V,rvect)*(h/R^4)) * y(1) - 2*(h/R^2)*y(4);
dydt(6) = -mu/(R^3) * y(3);

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