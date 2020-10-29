function  Rho = DensityModel(Phi, Rhomat, Rhofl)

% DENSITY MODEL implements the linear porosity-density relation
% INPUT Phi = Porosity
%       Rhomat = Density of the solid phase
%       Rhofl = Density of the fluid phase
% OUTUPT Rho = Density of saturated rock

% Written by Dario Grana (August 2020)

Rho = (1-Phi).*Rhomat+ Phi.*Rhofl;
