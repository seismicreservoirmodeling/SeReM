function C = SpatialCovariance2D(lmin, lmax, azim, theta, h, type)

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
        C = ExpCov(h, RadialCorrLength(lmin, lmax, azim, theta));    
    case 'gau'
        C = GauGov(h, RadialCorrLength(lmin, lmax, azim, theta));
    case 'sph'
        C = SphCov(h, RadialCorrLength(lmin, lmax, azim, theta));
end
        
