function [dstate] = coast_ODE(time, state, mu)

x = state(1);
y = state(2);
z = state(3);
dx = state(4);
dy = state(5);
dz = state(6);
r = norm([x y z]);
ddx = -mu*x/r^3;
ddy = -mu*y/r^3;
ddz = -mu*z/r^3;
dstate = [dx; dy; dz; ddx; ddy; ddz];

end