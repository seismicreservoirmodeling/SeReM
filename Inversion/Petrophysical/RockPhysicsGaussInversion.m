function [mupost, sigmapost, Ppost]  = RockPhysicsGaussInversion(mtrain, dtrain, mdomain, dcond, sigmaerr)

% ROCK PHYSICS GAUSSIAN INVERSION computes the posterior distribution of
% petrophysical properties conditioned on elastic properties assuming a
% Gaussian distribution.
% The joint distribution of the Bayesian inversion approach is estimated
% from a training dataset 
% INPUT mtrain = training dataset of petrophysical properties (ntrain, nm)
%       dtrain = training dataset of elastic properties (ntrain, nd)
%       mdomain = discretized domain of petrophysical properties 
%                 (generated using meshgrid)
%       dcond = measured data (nsamples, nd)
%       sigmaerr = covariance matrix of the error (nd, nd)
% OUTUPT mupost = posterior mean (nsamples x nv, 1)
%        sigmapost = posterior covariance matrix (nv, nv) 
%        Ppost = joint posterior distribution 

% initial parameters
nv = size(mtrain, 2);
ns = size(dcond,1);
datatrain = [mtrain dtrain];

% joint distribution
mjoint = mean(datatrain);
mum = mjoint(1:nv);
mud = mjoint(nv+1:end);
sjoint = cov(datatrain);
sm = sjoint(1:nv,1:nv);
sd = sjoint(nv+1:end,nv+1:end);
smd = sjoint(1:nv,nv+1:end);
sdm = sjoint(nv+1:end,1:nv);

% posterior distribution
mupost = zeros(ns,nv);
Ppost = zeros(ns,size(mdomain,1));
% posterior covariance matrix
sigmapost = sm-smd/(sd+sigmaerr)*sdm;
% [~,posdefcheck] = chol(sigmapost);
% [V,D]=eig(sigmapost);
% d=diag(D);
% d(d<=0)=eps;
% sigmapost= V*diag(d)*V';
% analytical solution
for i=1:ns
    % posterior mean
    mupost(i,:) = mum'+smd/(sd+sigmaerr)*(dcond(i,:)-mud)';
    % posterior PDF
    Ppost(i,:)=mvnpdf(mdomain, mupost(i,:), sigmapost)';
end
