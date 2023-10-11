function [xok, xvarok] = OrdinaryKriging(xcoord, dcoords, dvalues, xvar, l, type, angles)

% ORDINARY KRIGING computes the ordinary kriging estimate and variance 
% INPUT xcoord = coordinates of the location for the estimation (1, ndim)
%       dcoords = coordinates of the measurements (ns, ndim)
%       dvalues = values of the measurements (ns, 1)
%       xvar = prior variance
%       h = distance
%       l = correlation length, 1x1 for isotropic or 3x1 for anisotropic
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

% if isotropic
if length(l) == 1
    krigvect(1:nd,1) = xvar*SpatialCovariance1D(distvect,l,type);
    krigmatr(1:nd,1:nd) = xvar*SpatialCovariance1D(distmatr,l,type);
else % if anisotropic
    % Apply transformation in the coordinate system according to the variogram parameters
    coords = [xcoord; dcoords];
    ROT = rotx(angles(1)) * roty(angles(2)) * rotz(angles(3));
    coords = (ROT * coords')';       
    coords = coords./l;
    
    xdtemp = squareform(pdist(coords));
    distvect = xdtemp(2:end,1);
    distmatr = xdtemp(2:end,2:end);
    krigvect(1:nd,1) = xvar*SpatialCovariance1D(distvect,1,type);
    krigmatr(1:nd,1:nd) = xvar*SpatialCovariance1D(distmatr,1,type);

end

    
krigmatr(end,end) = 0;

% kriging weights
wkrig = krigmatr\krigvect;
    
% kriging mean
% xok = mean(dvalues)+sum(wkrig(1:end-1).*(dvalues-mean(dvalues)));
xok = sum(wkrig(1:end-1).*(dvalues));
% kriging variance
xvarok = xvar-sum(wkrig.*krigvect);