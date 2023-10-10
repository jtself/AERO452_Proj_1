function h = findh(r, mu, ecc, theta)
% given magnitude r, mu, eccentricity, and true anomaly, this function
% finds the magnitude of h, the specific angular momentum of ANY orbit

h = sqrt(r*mu*(1 + ecc* cos(theta)));