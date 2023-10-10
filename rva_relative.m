function [r_relx, v_relx, a_relx] = rva_relative(rA,vA,rB,vB)
%{
 This function is from Curtis Algorithm 7.1.
This uses the state vectors of s/c A (target) and s/c B (chaser) to find
the position, velocity, and acceration of B with respect to A (relative to
A) in the LVLH frame. 

% Justin Self 2023

Nomenclature:
INPUTS:
rA, vA = state vector of A (km, km/s)
rB, vB = state vector of B (km, km/s)

NOTES:
mu = grav parameter of planet (km3/s2)
hA = angular momentum vector of A (km2/s)
i, j, k = unit vectors along the x, y, and z axes of A's LVLH frame
QXx = Direction Cosine Matrix of the LVLH frame relative to geocentric
equatorial frame (GEF) (or ECI)
Omega = angular velocity of the LVLH frame (rad/s)
Omega_dot = angular acceleration of the LVLH frame (rad/s2)
aA, aB = absolute accelerations of A and B (km/s2)
r_rel = position of B relative to A in ECI
v_rel = velocity of B relative to A in ECI (km/s) 
a_rel = acceleration of B relative to A in ECI (km/s^2)

OUTPUTS: 
r_relx - position of B relative to A in the LVLH frame
v_relx = velocity of B relative to A in the LVLH frame
a_relx = acceleration of B relative to A in the LVLH frame

%}

mu = 398600; % grav param earth, km3/s2

% convert from ECI to LVLH
    % CALCULATE ANGULAR MOMENTUM VECTOR OF S/C A
    hA = cross(rA,vA); % km2/s

    % CALCULATE UNIT VECTORS OF LVLH FRAME
    i_hat = rA / norm(rA); 
    k_hat = hA / norm(hA); % mag of h never changes (assume no perturbations)
    j_hat = cross(k_hat,i_hat);
    
    % CALCULATE ORTHOGONAL DIRECTION COSINE MATRIX (QXx) using Eq. 7.11
    QXx = [i_hat'; j_hat'; k_hat'];
    
    % CALCULATE CAP_OMEGA and CAP_OMEGA_DOT FROM Eqns 7.5 and 7.6
    Omega = hA / (norm(rA)^2);
    Omega_dot = -2*dot(rA,vA) / norm(rA)^2 * Omega;

    % Calculate accelerations aA and aB
    aA = -mu * rA / norm(rA)^3;
    aB = -mu * rB / norm(rB)^3;

    % in ECI still
    
    % Find rho dot (r_rel) and rho doubledot (v_rel) (in ECI)
    r_rel = rB - rA;
    v_rel = vB - vA - cross(Omega,r_rel);

    % Calculate a_rel
    a_rel = aB - aA - cross(Omega_dot,r_rel) - cross(Omega,cross(Omega,r_rel))...
        - 2 * cross(Omega,v_rel);

    % Calculate r_relx, v_relx, and a_relx (all in LVLH frame)
    r_relx = QXx*r_rel;
    v_relx = QXx*v_rel;
    a_relx = QXx*a_rel;

end % function