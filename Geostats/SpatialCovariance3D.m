function C = SpatialCovariance3D(l, angles, theta, gamma, h, type)

% SPATIAL COVARIANCE 2D computes the 2D anisotropic spatial covariance 
% function 
% INPUT lmin = minimum correlation length
%       lmax = aaximum correlation length
%       azim = azimuth
%       theta = radial coordinate
%       h = distance
%       type = function ype ('exp', 'gau', 'sph')
% OUTPUT C = covariance

% Written by Dario Grana (August, 2020)

% covariance function
switch type
    case 'exp'
        C = ExpCov(h, RadialCorrLength(l, angles, theta, gamma));    
    case 'gau'
        C = GauCov(h, RadialCorrLength(l, angles, theta, gamma));
    case 'sph'
        C = SphCov(h, RadialCorrLength(l, angles, theta, gamma));
end
        
