function C = ExpCov(h,l)

% EXP COV computes the exponential covariance function
% INPUT h = distance
%       l = correlation length (or range)
% OUTPUT C = covariance

% Written by Dario Grana (August, 2020)

% covariance function
C = exp(-3*h/l);


