function [r_ECI,v_ECI] = r_and_v_from_COEs(RAAN,inc,w,h,Ecc,theta)
% Finds the r and v vectors in an Earth Centered Inertial reference frame
% Inputs
    % inc = inclination [rad]
    % w = argiment of perigee [rad]
    % RAAN = right ascension of the ascending node [rad]
    % h = angular momentum [km2/s]
    % Ecc = eccentricity
    % theta = true anomaly [rad]
% Outputs
    % r_vect = position vector in ECI frame [km]
    % v_vect = velocity vector in ECI frame [km/s]

% constant    
muEarth = 398600;

% Convert to rad
RAAN = deg2rad(RAAN);
inc = deg2rad(inc);
w = deg2rad(w);
theta = deg2rad(theta);

% determine perifocal to GEO rotation matrix
Cz_RAAN =        [  cos(RAAN)     sin(RAAN)     0;
                   -sin(RAAN)     cos(RAAN)     0;
                    0             0             1]; % 3

Cx_inc =         [  1             0             0;
                    0             cos(inc)      sin(inc);
                    0             -sin(inc)     cos(inc)]; % 1

Cz_w =           [  cos(w)        sin(w)        0;
                   -sin(w)        cos(w)        0;
                    0             0             1]; % 3

Qgeo2peri = Cz_w * Cx_inc * Cz_RAAN;
Qperi2geo = Qgeo2peri';

% find r and v (Perifocal)
r_PF = (((h^2)/muEarth)*(1/(1+Ecc*cos(theta))))*[cos(theta);sin(theta);0];
v_PF = (muEarth/h)*[-sin(theta); (Ecc+cos(theta)); 0];

% calc r and v relative to the geocentric reference frame
% GEO = ECI
r_ECI = Qperi2geo*r_PF;
v_ECI = Qperi2geo*v_PF;