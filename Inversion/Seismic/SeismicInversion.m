function [mmap, mlp, mup, Time] = SeismicInversion(Seis, TimeSeis, Vpprior, Vsprior, Rhoprior, sigmaprior, sigmaerr, wavelet, theta, nv)

% SEISMIC INVERSION computes the posterior distribution of elastic
% properties according to the Bayesian linearized AVO inversion (Buland and
% Omre, 2003)
% INPUT Seis = vector of seismic data of size (nsamples x nangles, 1)
%       TimeSeis = vector of seismic time of size (nsamples, 1)
%       Vpprior = vector of prior (low frequency) Vp model (nsamples+1, 1)
%       Vsprior = vector of prior (low frequency) Vs model (nsamples+1, 1)
%       Rhoprior = vector of prior (low frequency) density model (nsamples+1, 1)
%       sigmaprior = prior covariance matrix (nv*(nsamples+1),nv*(nsamples+1))
%       sigmaerr = covariance matrix of the error (nv*nsamples,nv*nsamples)
%       theta = vector of reflection angles (1,nangles)
%       nv = number of model variables
% OUTUPT mmap = MAP of posterior distribution (nv*(nsamples+1),1)
%        mlp = P2.5 of posterior distribution (nv*(nsamples+1),1)
%        mup = P97.5 of posterior distribution (nv*(nsamples+1),1)
%        Time = time vector of elastic properties (nsamples+1,1)

% Written by Dario Grana (August 2020)

% parameters
ntheta = length(theta);

% logarithm of the prior
logVp = log(Vpprior);
logVs = log(Vsprior);
logRho = log(Rhoprior);
mprior = [logVp; logVs; logRho];
nm = size(logVp,1);

% Aki Richards matrix
A = AkiRichardsCoefficientsMatrix(Vpprior, Vsprior, theta, nv);

% Differential matrix 
D = DifferentialMatrix(nm,nv);

% Wavelet matrix
W = WaveletMatrix(wavelet, nm, ntheta);

% forward operator
G = W*A*D;

% Bayesian Linearized AVO inverison analytical solution (Buland and Omre, 2003)
% mean of d
mdobs = G*mprior;
% covariance matrix
sigmadobs = G*sigmaprior*G'+sigmaerr;

% posterior mean
mpost = mprior+(G*sigmaprior)'*(sigmadobs\(Seis-mdobs));
% posterior covariance matrix
sigmapost = sigmaprior-(G*sigmaprior)'*(sigmadobs\(G*sigmaprior));

% statistical estimators posterior distribution
mmap = exp(mpost-diag(sigmapost));
mlp = exp(mpost-1.96*sqrt(diag(sigmapost)));
mup = exp(mpost+1.96*sqrt(diag(sigmapost)));

% time
dt = TimeSeis(2)-TimeSeis(1);
Time = (TimeSeis(1)-dt/2:dt:TimeSeis(end)+dt/2)';