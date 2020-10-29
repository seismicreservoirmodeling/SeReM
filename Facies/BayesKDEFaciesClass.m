function [fmap, fpost] = BayesKDEFaciesClass(data, dtrain, ftrain, fprior, domain)

% BAYES KDE FACIES CLASS computes the Bayesian facies classification of
% the data assuming a non-parametric distribution
% INPUT data = input data (ns,nv)
%       dtrain = training data (ntrain,nv)
%       ftrain = training facies (ntrain,1)
%       fprior = prior facies proportions (nf,1)
%       domain = discretized domain of input vairables 
%                (generated using meshgrid)
% OUTUPT facies maximum a posteriori (ns,1)
%        fpost = posterior facies probability (ns,nf)

% Written by Dario Grana (August 2020)

% initial parameters
ns = size(data,1);
nf = max(unique(ftrain));
nd = size(domain,1);

% joint distribution
Pjoint = zeros(nd,nf);
for k=1:nf
    Pjoint(:,k) = ksdensity(dtrain(ftrain==k,:),domain, 'Kernel', 'epanechnikov');
    Pjoint(:,k) = Pjoint(:,k)/sum(Pjoint(:,k));
end

% conditional distribution
fmap = zeros(ns,1);
fpost = zeros(ns,nf);
for i=1:ns
    [~,ind] = min(sum((domain-data(i,:)).^2,2));
    for k=1:nf
        fpost(i,k) = fprior(k)*Pjoint(ind,k);
    end
    % probability
    fpost(i,:) = fpost(i,:)/sum(fpost(i,:));
    % maximum a posteriori
    [~, fmap(i,:)] = max(fpost(i,:));
end