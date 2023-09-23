function [Seis, TimeSeis, Cpp] = SeismicModel1D (Vp, Vs, Rho, Time, theta, W, D)

% SeismicModel3D computes synthetic seismic data according to a linearized
% seismic model based on the convolution of a wavelet and the linearized
% approximation of Zoeppritz equations
% INPUT Vp = P-wave velocity profile
%       Vs = S-wave velocity profile
%       Rho = Density profile
%       theta = vector of reflection angles 
%       W = wavelet matrix
%       D = differential matrix
% OUTUPT Seis = vector of seismic data of size (nsamples x nangles, 1)
%        Time = seismic time  (nsamples, 1)
%        Cpp  = reflectivity coefficients matrix (nsamples x nangles, nsamples x nangles)

% Written by Dario Grana (August 2020)
% Modified by Dario Grana (June 2023)


% number of variables
nv = 3;

% logarithm of model variables
logVp = log(Vp);
logVs = log(Vs);
logRho = log(Rho);
m = [logVp; logVs; logRho];

% Aki Richards matrix
A = AkiRichardsCoefficientsMatrix(Vp, Vs, theta, nv);

% Reflectivity coefficients matrix
mder=D*m;
Cpp = A*mder;

% Seismic data matrix
Seis = W*Cpp;

% Time seismic measurements
TimeSeis = 1/2*(Time(1:end-1)+Time(2:end));

