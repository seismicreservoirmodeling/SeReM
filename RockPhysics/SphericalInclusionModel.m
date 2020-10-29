function [Vp, Vs, Rho, Ksat, Gsat] = SphericalInclusionModel(Phi, Rho, Kmat, Gmat, Kfl)

% SPHERICAL INCLUSION MODEL implements the inclusion model for spherical
% pores
% INPUT Phi = Porosity
%       Rho = Density of the saturated rock
%       Kmat = Bulk modulus of the solid phase
%       Gmat = Shear modulus of the solid phase
%       Kfl = Bulk modulus of the fluid phase
% OUTUPT Vp = P=wave velocity
%        Vs = S-wave velocity

% Written by Dario Grana (August 2020)

% elastic moduli of the dry rock
Kdry=4*Kmat*Gmat*(1-Phi)./(3*Kmat.*Phi+4*Gmat);
Gdry=Gmat*(9*Kmat+8*Gmat).*(1-Phi)./((9*Kmat+8*Gmat+6*(Kmat+2*Gmat).*Phi));

% Gassmann
[Ksat, Gsat] = GassmannModel(Phi, Kdry, Gdry, Kmat, Kfl);

% Velocities
Vp = sqrt((Ksat+4/3*Gsat)./Rho);
Vs = sqrt(Gsat./Rho);