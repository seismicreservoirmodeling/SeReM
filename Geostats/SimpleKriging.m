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
xdtemp = squareform(pdist([xcoord; dcoords]));
distvect = xdtemp(2:end,1);
distmatr = xdtemp(2:end,2:end);
krigvect = xvar*SpatialCovariance1D(distvect,l,type);
krigmatr = xvar*SpatialCovariance1D(distmatr,l,type);

% kriging weights
wkrig = krigmatr\krigvect;
    
% kriging mean
xsk = xmean+sum(wkrig.*(dvalues-xmean));
% kriging variance
xvarsk = xvar-sum(wkrig.*krigvect);