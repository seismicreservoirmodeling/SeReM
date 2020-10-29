function p = EpanechnikovKernel(x)

% EPANECHNIKOV KERNEL computes Epanechnikov Kernel 
% INPUT x = independent variable
% OUTUPT p = kernel

% Written by Dario Grana (August 2020)

p = (0.75*(1-x.^2)).*(x>=-1 & x<=1);