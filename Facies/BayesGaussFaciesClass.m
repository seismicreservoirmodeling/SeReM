function [fmap, fpost] = BayesGaussFaciesClass(data, fprior, muprior, sigmaprior)

% BAYES GAUSS FACIES CLASS computes the Bayesian facies classification of
% the data assuming a Gaussian distribution
% INPUT data = input data (ns,nv)
%       fprior = prior facies proportions (nf,1)
%       muprior = prior means on input variables (nf,nv)
%       sigmaprior = prior covariancce matrices on input variables (nv,nv,nf)
% OUTUPT fmap = facies maximum a posteriori (ns,1)
%        fpost = posterior facies probability (ns,nf)

% Written by Dario Grana (August 2020)

% initial parameters
ns = size(data,1);
nf = size(muprior,1);

% conditional probability
fmap = zeros(ns,1);
fpost = zeros(ns,nf);
for i=1:ns
    for k=1:nf
        fpost(i,k) = fprior(k)*mvnpdf(data(i,:),muprior(k,:),sigmaprior(:,:,k));
    end
    % probability
    fpost(i,:) = fpost(i,:)/sum(fpost(i,:));
    % maximum a posteriori
    [~, fmap(i,:)] = max(fpost(i,:));
end