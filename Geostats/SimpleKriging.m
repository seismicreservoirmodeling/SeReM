function [xsk, xvarsk] = SimpleKriging(xcoord, dcoords, dvalues, xmean, xvar, l, type)

% SIMPLE KRIGING computes the simple kriging estimate and variance 
% INPUT xcoord = coordinates of the location for the estimation (1, ndim)
%       dcoords = coordinates of the measurements (ns, ndim)
%       dvalues = values of the measurements (ns, 1)
%       xmean = prior mean
%       xvar = prior variance
%       h = distance
%       l = correlation length
%       type = function ype ('exp', 'gau', 'sph')
% OUTPUT xsk = kriging estimate
%        xvarsk = kriging variance

% Written by Dario Grana (August, 2020)

% kriging matrix and vector
nd = size(dcoords,1);

coords = [xcoord; dcoords];
xdtemp = squareform(pdist([xcoord; dcoords]));
distvect = xdtemp(2:end,1);
distmatr = xdtemp(2:end,2:end);

if length(l) == 1
    krigvect(1:nd,1) = xvar*SpatialCovariance1D(distvect,l,type);
    krigmatr(1:nd,1:nd) = xvar*SpatialCovariance1D(distmatr,l,type);
else
    angles_rad = atan2( coords(:,2).' - coords(:,2) , coords(:,1).' - coords(:,1) );
    angles_rad = triu(angles_rad) + triu(angles_rad)';
    angvect = angles_rad(2:end,1);
    angmatr = angles_rad(2:end,2:end);
    krigvect(1:nd,1) = xvar*SpatialCovariance2D(l(1), l(2), 0, angvect, distvect, type);
    krigmatr(1:nd,1:nd) = xvar*SpatialCovariance2D(l(1), l(2), 0, angmatr, distmatr, type);
end

% kriging weights
wkrig = krigmatr\krigvect;
    
% kriging mean
xsk = xmean+sum(wkrig.*(dvalues-xmean));
% kriging variance
xvarsk = xvar-sum(wkrig.*krigvect);