function [xok, xvarok] = OrdinaryKriging(xcoord, dcoords, dvalues, xvar, l, type)

% ORDINARY KRIGING computes the ordinary kriging estimate and variance 
% INPUT xcoord = coordinates of the location for the estimation (1, ndim)
%       dcoords = coordinates of the measurements (ns, ndim)
%       dvalues = values of the measurements (ns, 1)
%       xvar = prior variance
%       h = distance
%       l = correlation length
%       type = function ype ('exp', 'gau', 'sph')
% OUTPUT xok = kriging estimate
%        xvarok = kriging variance

% Written by Dario Grana (August, 2020)

% kriging matrix and vector
nd = size(dcoords,1);
krigmatr = ones(nd+1,nd+1);
krigvect = ones(nd+1,1);
xdtemp = squareform(pdist([xcoord; dcoords]));
distvect = xdtemp(2:end,1);
distmatr = xdtemp(2:end,2:end);
krigvect(1:nd,1) = xvar*SpatialCovariance1D(distvect,l,type);
krigmatr(1:nd,1:nd) = xvar*SpatialCovariance1D(distmatr,l,type);
krigmatr(end,end) = 0;

% kriging weights
wkrig = krigmatr\krigvect;
    
% kriging mean
% xok = mean(dvalues)+sum(wkrig(1:end-1).*(dvalues-mean(dvalues)));
xok = sum(wkrig(1:end-1).*(dvalues));
% kriging variance
xvarok = xvar-sum(wkrig.*krigvect);