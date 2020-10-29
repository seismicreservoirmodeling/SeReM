function  [Ksat, Gsat] = GassmannModel(Phi, Kdry, Gdry, Kmat, Kfl)

% GASSMANN MODEL implements Gassmann's equations 
% INPUT Phi = Porosity
%       Kdry = Bulk modulus of dry rock
%       Gdry = Shear modulus of dry rock
%       Kmat = Bulk modulus of solid phase
%       Kfl = Bulk modulus of fluid rock
% OUTUPT Ksat = Bulk modulus of saturated rock
%        Gsat = Shear modulus of saturated rock

% Written by Dario Grana (August 2020)

% Bulk modulus of saturated rock
Ksat = Kdry + ((1-Kdry./Kmat).^2)./(Phi./Kfl+(1-Phi)./Kmat-Kdry./(Kmat.^2));
% Shear modulus of saturated rock
Gsat = Gdry;
