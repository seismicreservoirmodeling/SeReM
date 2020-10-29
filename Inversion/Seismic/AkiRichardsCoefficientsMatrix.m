function A = AkiRichardsCoefficientsMatrix(Vp, Vs, theta, nv)

% AKI RICHARDS COEFFICIENTS MATRIX computes the Aki Richards coefficient
% matrix
% INPUT Vp = P-wave velocity profile
%       Vs = S-wave velocity profile
%       theta = vector of reflection angles
%       nv = number of model variables 
% OUTUPT A = Aki Richards coefficients matrix

% Written by Dario Grana (August 2020)

% initial parameters
nsamples = length(Vp);
ntheta = length(theta);
A = zeros((nsamples-1)*ntheta, nv*(nsamples-1));

% average velocities at the interfaces
avgVp = 1/2*(Vp(1:end-1)+Vp(2:end));
avgVs = 1/2*(Vs(1:end-1)+Vs(2:end));

% reflection coefficients (Aki Richards linearized approximation)
for i=1:ntheta
    cp = 1/2*(1+tand(theta(i)).^2)*ones(nsamples-1,1);
    cs = -4*avgVs.^2./avgVp.^2*sind(theta(i))^2;
    cr = 1/2*(1-4*avgVs.^2./avgVp.^2*sind(theta(i))^2);
    Acp = diag(cp);
    Acs = diag(cs);
    Acr = diag(cr);   
    A((i-1)*(nsamples-1)+1:i*(nsamples-1),:) = [Acp, Acs, Acr];
end
