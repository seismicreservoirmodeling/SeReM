function C = GauCov(h,l)

% GAU COV computes the Gaussian covariance function
% INPUT h = distance
%       l = correlation length (or range)
% OUTPUT C = covariance

% Written by Dario Grana (August, 2020)

% covariance function
C = exp(-3*h.^2/l.^2);

