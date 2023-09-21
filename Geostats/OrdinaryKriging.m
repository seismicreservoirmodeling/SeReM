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

coords = [xcoord; dcoords];
xdtemp = squareform(pdist(coords));
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

    
krigmatr(end,end) = 0;

% kriging weights
wkrig = krigmatr\krigvect;
    
% kriging mean
% xok = mean(dvalues)+sum(wkrig(1:end-1).*(dvalues-mean(dvalues)));
xok = sum(wkrig(1:end-1).*(dvalues));
% kriging variance
xvarok = xvar-sum(wkrig.*krigvect);