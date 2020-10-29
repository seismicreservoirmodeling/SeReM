function  [Vp, Vs] = StiffsandModel(Phi, Rho, Kmat, Gmat, Kfl, critporo, coordnum, press)

% STIFF SAND MODEL implements Dvorkin's soft sand model
% INPUT Phi = Porosity
%       Rho = Density of the saturated rock
%       Kmat = Bulk modulus of the solid phase
%       Gmat = Shear modulus of the solid phase
%       Kfl = Bulk modulus of the fluid phase
%       critporo = critical porosity
%       coordnum = coordination number
%       press = effective pressure in GPA
% OUTUPT Vp = P=wave velocity
%        Vs = S-wave velocity

% Written by Dario Grana (August 2020)

% Hertz-Mindlin
Poisson = (3*Kmat-2*Gmat)./(6*Kmat+2*Gmat);
KHM = ((coordnum^2*(1-critporo)^2*Gmat.^2*press)./(18*pi^2*(1-Poisson).^2)).^(1/3);
GHM = (5-4*Poisson)./(10-5*Poisson).*((3*coordnum^2*(1-critporo)^2*Gmat.^2*press)./(2*pi^2*(1-Poisson).^2)).^(1/3);

% Modified Hashin-Shtrikmann upper bounds
Kdry = 1./((Phi/critporo)./(KHM+4/3*Gmat)+(1-Phi/critporo)./(Kmat+4/3*Gmat))-4/3*Gmat;
psi = (9*Kmat+8*Gmat)./(Kmat+2*Gmat);
Gdry = 1./((Phi/critporo)./(GHM+1/6*psi.*Gmat)+(1-Phi/critporo)./(Gmat+1/6*psi.*Gmat))-1/6*psi.*Gmat;

% Gassmann
[Ksat, Gsat] = GassmannModel(Phi, Kdry, Gdry, Kmat, Kfl);

% Velocities
Vp = sqrt((Ksat+4/3*Gsat)./Rho);
Vs = sqrt(Gsat./Rho);
