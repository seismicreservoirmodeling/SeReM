function  [Vp, Vs] = SoftsandModel(Phi, Rho, Kmat, Gmat, Kfl, critporo, coordnum, press)

% SOFT SAND MODEL implements Dvorkin's soft sand model
% INPUT Phi = Porosity
%       Rho = Density of the saturated rock
%       Kmat = Bulk modulus of the solid phase
%       Gmat = Shear modulus of the solid phase
%       Kfl = Bulk modulus of the fluid phase
%       critporo = critical porosity
%       coordnum = coordination number
%       pressure = effective pressure in GPA
% OUTUPT Vp = P=wave velocity
%        Vs = S-wave velocity

% Written by Dario Grana (August 2020)

% Hertz-Mindlin
Poisson = (3*Kmat-2*Gmat)./(6*Kmat+2*Gmat);
KHM = ((coordnum^2*(1-critporo)^2*Gmat.^2*press)./(18*pi^2*(1-Poisson).^2)).^(1/3);
GHM = (5-4*Poisson)./(10-5*Poisson).*((3*coordnum^2*(1-critporo)^2*Gmat.^2*press)./(2*pi^2*(1-Poisson).^2)).^(1/3);
% f = friction
% GHM = (2+3*f-Poisson*(1+3f))./(10-5*Poisson).*((3*coordnumber^2*(1-criticalporo)^2*Gmat.^2*pressure)./(2*pi^2*(1-Poisson).^2)).^(1/3);

% Modified Hashin-Shtrikmann lower bounds
Kdry = 1./((Phi/critporo)./(KHM+4/3*GHM)+(1-Phi/critporo)./(Kmat+4/3*GHM))-4/3*GHM;
psi = (9*KHM+8*GHM)./(KHM+2*GHM);
Gdry = 1./((Phi/critporo)./(GHM+1/6*psi.*GHM)+(1-Phi/critporo)./(Gmat+1/6*psi.*GHM))-1/6*psi.*GHM;

% Gassmann
[Ksat, Gsat] = GassmannModel(Phi, Kdry, Gdry, Kmat, Kfl);

% Velocities
Vp = sqrt((Ksat+4/3*Gsat)./Rho);
Vs = sqrt(Gsat./Rho);
