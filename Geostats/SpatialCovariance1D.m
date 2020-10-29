function C = SpatialCovariance1D(h, l, type)

% SPATIAL COVARIANCE 1D computes the 1D spatial covariance function 
% INPUT l = correlation length
%       h = distance
%       type = function ype ('exp', 'gau', 'sph')
% OUTPUT C = covariance

% Written by Dario Grana (August, 2020)

% covariance function
switch type
    case 'exp'
        C = ExpCov(h,l);    
    case 'gau'
        C = GauGov(h,l);
    case 'sph'
        C = SphCov(h,l);
end
        
