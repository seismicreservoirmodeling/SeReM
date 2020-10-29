function msim = CorrelatedSimulation(mprior, sigma0, sigmaspace)

% CORRELATED SIMULATION generates 1D stochastic realizations of correlated
% multiple random variables with a spatial correlation model
% INPUT mprior = prior trend (nsamples, nvariables)
%       sigma0 = stationary covariance matrix (nvariables, nvariables)
%       sigmaspace =spatial covariance matrix (nsamples, nsamples)
% OUTPUT msim = stochastic realization (nsamples, nvariables)

% Written by Dario Grana (August, 2020)

% initial parameters
nm = size(mprior,2);
ns = size(mprior,1);

% spatial covariance matrix
sigma = kron(sigma0, sigmaspace);

% spatially correlated realization
mreal = mvnrnd(mprior(:), sigma)';
msim = zeros(ns,nm);
for i=1:nm
    msim(:,i)=mreal((i-1)*ns+1:i*ns);
end
 