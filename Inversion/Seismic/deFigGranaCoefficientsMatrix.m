function A = deFigGranaCoefficientsMatrix(Vp, Vs, theta, nv)
% deFigGranaCoefficientsMatrix computes the de Figueiredo and Grana coefficient
% matrix
% INPUT Vp = P-wave velocity profile
%       Vs = S-wave velocity profile
%       theta = vector of reflection angles
%       nv = number of model variables 
% OUTUPT A = Aki Richards coefficients matrix

% Written by Dario Grana (August 2020)
% Edited with sparse by Jorlivan Correa (December 2022)

% initial parameters
nsamples = length(Vp);
ntheta = length(theta);
A = sparse((nsamples-1)*ntheta, nv*(nsamples-1));

% average velocities at the interfaces
avgVp = 1/2*(Vp(1:end-1)+Vp(2:end));
avgVs = 1/2*(Vs(1:end-1)+Vs(2:end));

% reflection coefficients (Aki Richards linearized approximation)
for i=1:ntheta
    cp = (0.5 - (18/5)*(avgVs.^2./avgVp.^2)*sind(theta(i))^2 + (2/5)*tand(theta(i))^2 ).*ones(nsamples-1,1);
    cps = +4*avgVs.^2./avgVp.^2*sind(theta(i))^2;   
    Acp = diag(sparse(cp));
    Acps = diag(sparse(cps));
    A((i-1)*(nsamples-1)+1:i*(nsamples-1),:) = [Acp, Acps];
end