function [Seis, TimeSeis] = SeismicModel (Vp, Vs, Rho, Time, theta, wavelet)

% SEISMIC MODEL computes synthetic seismic data according to a linearized
% seismic model based on the convolution of a wavelet and the linearized
% approximation of Zoeppritz equations
% INPUT Vp = P-wave velocity profile
%       Vs = S-wave velocity profile
%       Rho = Density profile
%       theta = vector of reflection angles 
%       wavelet = wavelet
% OUTUPT Seis = vector of seismic data of size (nsamples x nangles, 1)
%        Time = seismic time  (nsamples, 1)

% Written by Dario Grana (August 2020)

% initial parameters
ntheta = length(theta);
nm = length(Vp);

% number of variables
nv = 3;

% logarithm of model variables
logVp = log(Vp);
logVs = log(Vs);
logRho = log(Rho);
m = [logVp; logVs; logRho];

% Aki Richards matrix
A = AkiRichardsCoefficientsMatrix(Vp, Vs, theta, nv);

% Differential matrix 
D = DifferentialMatrix(nm,nv);
mder=D*m;

% Reflectivity coefficients matrix
Cpp = A*mder;

% Wavelet matrix
W = WaveletMatrix(wavelet, nm, ntheta);

% Seismic data matrix
Seis = W*Cpp;

% Time seismic measurements
TimeSeis = 1/2*(Time(1:end-1)+Time(2:end));

