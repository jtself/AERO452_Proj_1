function [dv0_PLUS_start_burn,dvf_MINUS_off_burn,deltaV,deltaV_afterBurn] = cw_twoimpulse(dr,drf,dv0,period,t)

% This function uses the equation form of the Clohessy-Wiltshire equations
% to calculate position and velocity (relative) with known input values.
% CIRCULAR ORBIS ONLY.

% INPUT time, period (in seconds), and delta_r vector (3x1)
% OUTPUT: delta_v +

n = 2*pi / period; % mean motion

% Matrix Method
phiRR = [4 - 3*cos(n*t),         0,     0;
         6*(sin(n*t) - (n*t)),   1,     0;
         0,                      0,     cos(n*t)];

phiRV = [(1/n)*sin(n*t),         (2/n)*(1-cos(n*t)) ,           0;
         (2/n)*(cos(n*t) - 1),   (1/n)*(4*sin(n*t)-(3*n*t)),    0;
         0,                      0,                             (1/n)*sin(n*t)];

phiVR = [3*n*sin(n*t),         0,           0;
         6*n*(cos(n*t) - 1),   0,           0;
         0,                    0,           -n*sin(n*t)];

phiVV = [cos(n*t),          2*sin(n*t) ,            0;
         -2*sin(n*t),       4*cos(n*t)-3,           0;
         0,                      0,                 cos(n*t)];



% Final outputs
dv0_PLUS_start_burn = inv(phiRV)*(drf + (-phiRR*dr));

dvf_MINUS_off_burn = phiVR*dr + phiVV*dv0_PLUS_start_burn;

deltaV = norm(dv0_PLUS_start_burn - dv0) + norm(dvf_MINUS_off_burn);
deltaV_afterBurn = dv0_PLUS_start_burn + dv0;

end % function