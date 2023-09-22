function [xok, xvarok] = OrdinaryKriging(xcoord, dcoords, dvalues, xvar, l, type)

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
else
    % if anisotropic
    azim_rad = atan2( coords(:,2).' - coords(:,2) , coords(:,1).' - coords(:,1) );
    azim_rad = triu(azim_rad) + triu(azim_rad)';
    azimvect = azim_rad(2:end,1);
    azimmatr = azim_rad(2:end,2:end);
    
    if length(l) ==3
        r = sqrt( coords(:,1).^2 + coords(:,2).^2 );
        dip_rad = atan2( r.' - r , coords(:,3).' - coords(:,3) );
        dip_rad = triu(dip_rad) + triu(dip_rad)';
        dipvect = dip_rad(2:end,1);
        dipmatr = dip_rad(2:end,2:end);
    else
        dipvect = pi/2*ones(size(azimvect));
        dipmatr = pi/2*ones(size(azimmatr));
    end
    
    krigvect(1:nd,1) = xvar*SpatialCovariance3D(l, [0,pi/2], azimvect, dipvect, distvect, type);
    krigmatr(1:nd,1:nd) = xvar*SpatialCovariance3D(l, [0,pi/2], azimmatr, dipmatr, distmatr, type);
end

    
krigmatr(end,end) = 0;

% kriging weights
wkrig = krigmatr\krigvect;
    
% kriging mean
% xok = mean(dvalues)+sum(wkrig(1:end-1).*(dvalues-mean(dvalues)));
xok = sum(wkrig(1:end-1).*(dvalues));
% kriging variance
xvarok = xvar-sum(wkrig.*krigvect);