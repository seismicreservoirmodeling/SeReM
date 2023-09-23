function [mmap, mtrans, strans, Time] = SeismicInversion3D(Seis, TimeSeis, Vpprior, Vsprior, Rhoprior, sigmaprior, sigmaerr, wavelet, theta, nv)

% SeismicInversion3D computes the posterior distribution of elastic
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
%        mtrans = transformed mean (logGaussian to Gaussian) (nv*(nsamples+1),1)
%        strans = transformed variance (logGaussian to Gaussian) (nv*(nsamples+1),1)
%        Time = time vector of elastic properties (nsamples+1,1)

% Written by Dario Grana (August 2020)
% Modified by Dario Grana (June 2023)

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

% posterior distribution
mtrans =  exp(mpost+1/2*diag(sigmapost));
% strans = sigmapost;
% for i=1:size(sigmapost,1)
%     for j=1:size(sigmapost,2)
%         strans(i,j) = exp(mpost(i)+mpost(j)+1/2*(sigmapost(i,i)+sigmapost(j,j))).*(exp(sigmapost(i,j))-1);
%     end
% end
strans = exp(2*mpost+diag(sigmapost)).*(exp(diag(sigmapost))-1);

% time
dt = TimeSeis(2)-TimeSeis(1);
Time = (TimeSeis(1)-dt/2:dt:TimeSeis(end)+dt)';