function l = RadialCorrLength(lmin, lmax, azim, theta)

% RADIAL CORR LENGTH computes the radial correlation length 
% INPUT lmin = minimum correlation length
%       lmax = aaximum correlation length
%       azim = azimuth
%       theta = radial coordinate
% OUTPUT l = radial correlation length

% Written by Dario Grana (August, 2020)

% covariance function
l = sqrt((lmin^2*lmax^2)./(lmax^2*(sin(azim-theta)).^2+lmin^2*(cos(azim-theta)).^2));
