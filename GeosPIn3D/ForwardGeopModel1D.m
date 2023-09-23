function SeisPred = ForwardGeopModel1D(Phi, Clay, Sw, regcoef, Time, theta, wavelet, nm, nv, nsim)

% ForwardGeopModel3D computes predicted seismic data
% INPUT Phi = porosity model (nxl x nil x nm)
%       Clay = clay model (nxl x nil x nm)
%       Sw = saturation model (nxl x nil x nm)
%       regcoef = Regression coefficients rock physics model (3x4)
%       Time = time vector (nd x 1)
%       theta = vector of reflection angles 
%       wavelet = wavelet vector 
%       nm = number of samples
%       nv = number of model variables
%       nsim = number of simulations
% OUTPUT SeisPred = Predicted seismic (nxl x nil x nd)

% Written by Dario Grana (June 2023)

nd = nm-1;
ntheta = length(theta);
SeisPred = zeros(ntheta*nd, nsim);
[Vp, Vs, Rho] = LinearizedRockPhysicsModel(Phi, Clay, Sw, regcoef);
D = DifferentialMatrix(nm,nv);
W = WaveletMatrix(wavelet, nm, ntheta);
for k=1:nsim
    [SeisPred(:,k), ~] = SeismicModel1D (Vp(:,k), Vs(:,k), Rho(:,k), Time, theta, W, D);
end
