function dydt = VBar_approach_CW(time,y,n,vc,mu)

dydt(1) = 0;
dydt(2) = vc;
dydt(3) = 0;

dydt(4) = (-n*2*vc);
dydt(5) = 0;
dydt(6) = 0;


% CHIEF, target
% Just propogating forward like normal using acceleration function (pig function).

R = norm([y(7) y(8) y(9)]);

dydt(7) = y(10);
dydt(8) = y(11);
dydt(9) = y(12);
dydt(10) = -mu * y(7) / R^3;
dydt(11) = -mu * y(8) / R^3;
dydt(12) = -mu * y(9) / R^3;

dydt = dydt';

end