function Ppost = RockPhysicsKDEInversion(mtrain, dtrain, mdomain, ddomain, dcond, hm, hd)

% ROCK PHYSICS KDE INVERSION computes the posterior distribution of
% petrophysical properties conditioned on elastic properties assuming a
% non-parametric distribution.
% The joint distribution of the Bayesian inversion approach is estimated
% from a training dataset using Kernel Density Estimation
% INPUT mtrain = training dataset of petrophysical properties (ntrain, nm)
%       dtrain = training dataset of elastic properties (ntrain, nd)
%       mdomain = discretized domain of petrophysical properties (ndiscr, nm)
%       ddomain = discretized domain of elastic properties (ndiscr, nd)
%       dcond = measured data (nsamples, nd)
%       hm = kernel bandwidths  hs of petrophysical properties (nm, 1)
%       hd = kernel bandwidths  of elastic properties (nd, 1)
% OUTUPT Ppost = joint posterior distribution 

% Written by Dario Grana (August 2020)

% number of training datapoint
nt = size(mtrain,1);
% number of data points
ns = size(dcond,1);

% number of model variables
nm = size(mtrain,2);
% number of data variables
nd = size(dtrain,2);

% multidimentional grids for model variables
vecm = num2cell(mdomain,1);
nem = numel(vecm);
mgrid = cell(1,nem);
[mgrid{1:end}] = ndgrid(vecm{1:end});
mgrid = cat(nem+1, mgrid{:});
mgrid = reshape(mgrid,[],nem);
% multidimentional grids for data variables
vecd = num2cell(ddomain,1);
ned = numel(vecd);
dgrid = cell(1,ned);
[dgrid{1:end}] = ndgrid(vecd{1:end});
dgrid = cat(ned+1, dgrid{:});
dgrid = reshape(dgrid,[],ned);

% kernel density estimation
mgridrep = zeros(nt,nm,size(mgrid,1));
for i=1:nt
    mgridrep(i,:,:) = mgrid';
end
datamrep = repmat(mtrain,[1,1,size(mgrid,1)]);
hmmat = repmat(hm,[nt,1,size(mgrid,1)]);
prodm = squeeze(prod(EpanechnikovKernel((mgridrep-datamrep)./hmmat),2));
dgridrep = zeros(nt,nd,size(dgrid,1));
for i=1:nt
    dgridrep(i,:,:) = dgrid';
end
datadrep = repmat(dtrain,[1,1,size(dgrid,1)]);
hdmat = repmat(hd,[nt,1,size(dgrid,1)]);
prodd = squeeze(prod(EpanechnikovKernel((dgridrep-datadrep)./hdmat),2));

% joint distribution
Pjoint = prodm'*prodd/(nt*prod(hm)*prod(hd));
Pjoint=Pjoint/sum(Pjoint(:));

% posterior distribution
Ppost=zeros(ns,size(mgrid,1));
for i=1:ns
    [~,indcond]=min(sum((dgrid-repmat(dcond(i,:),size(dgrid,1),1)).^2, 2));
    if sum(Pjoint(:,indcond))>0
        Ppost(i,:)=Pjoint(:,indcond)./sum(Pjoint(:,indcond));
    else
        Ppost(i,:)=Pjoint(:,indcond);
    end
end


