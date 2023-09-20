function l = RadialCorrLength(l1, l2, l3, angles, theta)
%function l = RadialCorrLength(lmin, lmax, azim, theta)

% RADIAL CORR LENGTH computes the radial correlation length 
% INPUT lmin = minimum correlation length
%       lmax = aaximum correlation length
%       azim = azimuth
%       theta = radial coordinate
% OUTPUT l = radial correlation length

% Written by Dario Grana (August, 2020)

% covariance function
l = sqrt( (lmin^2*lmax^2) ./ (lmax^2*(sin(azim-theta)).^2 + lmin^2*(cos(azim-theta)).^2) );

%azim = angles(1);
%dip = angles(2);
%l = sqrt( (l1^2*l2^2*l3^2) ./ ( l3^2*(l2^2*(cos(azim-theta)).^2 + l1^2*(sin(azim-theta)).^2)*(cos(dip))^2 + l1^2*l2^2*sin(dip)^2 ) );
