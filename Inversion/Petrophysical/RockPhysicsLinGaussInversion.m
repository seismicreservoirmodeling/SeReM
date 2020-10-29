function [mupost, sigmapost, Ppost]  = RockPhysicsLinGaussInversion(mum, sm, G, mdomain, dcond, sigmaerr)

% ROCK PHYSICS LINEAR GAUSSIAN INVERSION computes the posterior 
% distribution petrophysical properties conditioned on elastic properties 
% assuming a Gaussian distribution and a linear rock physics model
% INPUT mum = prior mean of petrophysical propertiees (nv, 1)
%       sm = prior covariance matrix of petrophysical propertiees (nv, nv)
%       G = rock physics operator matrix
%       mdomain = discretized domain of petrophysical properties 
%                 (generated using meshgrid)
%       dcond = measured data (nsamples, nd)
%       sigmaerr = covariance matrix of the error (nd, nd)
% OUTUPT mupost = posterior mean (nsamples x nv, 1)
%        sigmapost = posterior covariance matrix (nv, nv) 
%        Ppost = joint posterior distribution 

% Written by Dario Grana (August 2020)

% initial parameters
nv = size(mum,2);
ns = size(dcond,1);

% analytical calculations
mud = G*mum';
sd = G*sm*G';
smd = sm*G';
sdm = G*sm;

% posterior distribution
mupost = zeros(ns,nv);
Ppost = zeros(ns,size(mdomain,1));
% posterior covariance matrix
sigmapost = sm-smd/(sd+sigmaerr)*sdm;
% analytical solution
for i=1:ns
    % posterior mean
    mupost(i,:) = mum'+smd/(sd+sigmaerr)*(dcond(i,:)-mud')';
    % posterior PDF
    Ppost(i,:)=mvnpdf(mdomain, mupost(i,:), sigmapost)';
end
