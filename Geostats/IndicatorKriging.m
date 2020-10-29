function [ikp, ikmap] = IndicatorKriging(xcoord, dcoords, dvalues, nf, pprior, l, type)

% INDICATOR KRIGING computes the indicator kriging estimate and variance 
% INPUT xcoord = coordinates of the location for the estimation (1, ndim)
%       dcoords = coordinates of the measurements (ns, ndim)
%       dvalues = values of the measurements (ns, 1)
%       nf = number of possible outcomes (e.g. number of facies)
%       pprior = prior probability (1,nf)
%       h = distance 
%       l = correlation length
%       type = function ype ('exp', 'gau', 'sph')
% OUTPUT ikp = indicator kriging probability
%        ikmap = maximum a posteriori of indicator kriging probability

% Written by Dario Grana (August, 2020)

% indicator variables
nd = size(dvalues,1);
indvar = zeros(nd,nf);
for i=1:nd
    indvar(i,dvalues(i)) = 1;
end

% kriging weights
xdtemp = squareform(pdist([xcoord; dcoords]));
distvect = xdtemp(2:end,1);
distmatr = xdtemp(2:end,2:end);
varprior = zeros(1,nf);
krigvect = zeros(nd,nf);
krigmatr = zeros(nd,nd,nf);
wkrig = zeros(nd,nf);
for j=1:nf
    varprior(:,j) = pprior(j).*(1-pprior(j));
    krigvect(:,j) = varprior(:,j)*SpatialCovariance1D(distvect,l,type);
    krigmatr(:,:,j) = varprior(:,j)*SpatialCovariance1D(distmatr,l,type);
    wkrig(:,j) = krigmatr(:,:,j)\krigvect(:,j);
end

% indicator kriging probability
ikp = zeros(1,nf);
for j=1:nf
   ikp(1,j) = pprior(j)+sum(wkrig(:,j).*(indvar(:,j)-pprior(j)));
end
[~, ikmap] = max(ikp);