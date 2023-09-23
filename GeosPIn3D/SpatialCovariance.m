function sigma = SpatialCovariance(corrlength, dt, nm, sigma0)

% SpatialCovariance computes the spatial covaraince matrix
% INPUT corrlength = correlation length
%       dt = time interval 
%       nm = numbr of samples
%       sigma0 = covariance matrix
% OUTUPT sigmaprior = spatial covaraince matrix

% Written by Dario Grana (August 2020)

trow = repmat(0:dt:(nm-1)*dt,nm,1);
tcol = repmat((0:dt:(nm-1)*dt)',1,nm);
tdis = abs(trow-tcol);
sigmatime = exp(-(tdis./corrlength).^2);
sigma = kron(sigma0, sigmatime);
