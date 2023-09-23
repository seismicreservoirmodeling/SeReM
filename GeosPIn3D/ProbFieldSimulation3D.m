function [Phisim, Claysim, Swsim] = ProbFieldSimulation3D (vertcorr,horcorr,phiprior,clayprior, swprior, stdpetro, corrpetro, nxl, nil, nm, nsim)

% ProbFieldSimulation3D simulates spatially correlated realizations of 
% petrophysical properties using PFS simulations
% INPUT vertcorr = vertical correlation parameter
%       horcorr = horizontal correlation parameter
%       phprior = prior porosity mean
%       clayprior = prior clay mean
%       swprior = prior saturation mean
%       stdpetro = vector with prior standard deviation of petrophysical 
%                   properties (3 x 1)
%       corrpetro = prior correlation matrix of petrophysical properties (3x3)
%       nxl = number of crossline
%       nil - number of inline
%       nm = number of samples
%       nsim = number of simulations
% OUTPUT Phisim = Simulated porosity (nxl x nil x nm x nsim)
%       Claysim = Simulated clay (nxl x nil x nm x nsim)
%       Swsim = Simulated saturation (nxl x nil x nm x nsim)

% Written by Leandro de Figueiredo (March 2018)
% Modified by Dario Grana (June 2023)

Phisim = zeros(nxl,nil,nm,nsim);
Claysim = zeros(nxl,nil,nm,nsim);
Swsim = zeros(nxl,nil,nm,nsim);
corrfun = CorrelationFunction3D(vertcorr,horcorr,nxl, nil, nm);
% corrpetro = (diag(sqrt(diag(sigmapetro)))\sigmapetro)/diag(sqrt(diag(sigmapetro)));
for k=1:nsim
    uncorrsim = mvnrnd([0 0 0], corrpetro, nxl*nil*nm);
    noisephi = reshape(uncorrsim(:,1), nxl, nil, nm);
    noiseclay = reshape(uncorrsim(:,2), nxl, nil, nm);
    noisesw = reshape(uncorrsim(:,3), nxl, nil, nm);
    Phisim(:,:,:,k) = phiprior + stdpetro(1) *real(ifftn(sqrt(abs(fftn(corrfun))).*fftn(noisephi)));
    Claysim(:,:,:,k) = clayprior + stdpetro(2) *real(ifftn(sqrt(abs(fftn(corrfun))).*fftn(noiseclay)));
    Swsim(:,:,:,k) = swprior + stdpetro(3) *real(ifftn(sqrt(abs(fftn(corrfun))).*fftn(noisesw)));
end
