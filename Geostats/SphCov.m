function C = SphCov(h,l)

% SPH COV computes the spherical covariance function
% INPUT h = distance
%       l = correlation length (or range)
% OUTPUT C = covariance

% Written by Dario Grana (August, 2020)

% covariance function
C=zeros(size(h));
C(h<=l) = 1-3/2*h(h<=l)/l+1/2*h(h<=l).^3/l^3;

