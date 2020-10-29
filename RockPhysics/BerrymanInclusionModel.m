function [Vp, Vs] = BerrymanInclusionModel(Phi, Rho, Kmat, Gmat, Kfl, Ar)

% BERRYMAN INCLUSION MODEL implements Berryman's inclusion model for
% prolate and oblate spheroids
% INPUT Phi = Porosity
%       Rho = Density of the saturated rock
%       Kmat = Bulk modulus of the solid phase
%       Gmat = Shear modulus of the solid phase
%       Kfl = Bulk modulus of the fluid phase
%       Ar = Aspect ratio
% OUTUPT Vp = P=wave velocity
%        Vs = S-wave velocity

% Written by Dario Grana (August 2020)

% inclusion properties 
Kinc=Kfl;
Ginc = 0;

% Berryman's formulation
Poisson= (3*Kmat-2*Gmat)./(2*(3*Kmat+Gmat));
theta = Ar./(1-Ar.^2).^(3./2) .* (acos(Ar)-Ar.*sqrt(1-Ar.^2)); 
g = Ar.^2 ./(1-Ar.^2) .* (3.*theta-2);
R=(1-2.*Poisson)./(2-2.*Poisson);
A = (Ginc./Gmat)-1;
B = 1./3.*(Kinc./Kmat - Ginc./Gmat);
F1 = 1 + A .* (3./2.*(g + theta) - R .* (3./2.*g + 5./2.*theta - 4./3));
F2 = 1 + A .* (1 + 3./2.*(g + theta) - R./2.*(3.*g + 5 .* theta)) + B.*(3 - 4.*R)+A./2.*(A + 3.*B).*(3 - 4.*R).*(g + theta - R.*(g - theta + 2.*theta.^2));
F3 = 1 + A.*(1 - (g + 3/2 .* theta) + R.*(g + theta));
F4 = 1 + A./4.*(g + 3.*theta - R.*(g - theta));
F5 = A.*(R.*(g + theta - 4./3)-g) + B.*theta.*(3 - 4 .* R);
F6 = 1 + A.*(1 + g - R.*(theta+g)) + B.*(1 - theta) .* (3 - 4 .* R);
F7 = 2 + A./4 .* (9.*theta +3.*g- R.*(5.*theta + 3.*g)) + B.*theta.*(3 - 4.*R);
F8 = A.*(1 - 2.*R + g./2.*(R - 1) + theta./2.*(5.*R - 3)) + B.*(1 - theta).*(3 - 4.*R);
F9 = A.*(g.*(R - 1) - R.*theta) + B.*theta.*(3 - 4.*R); 
Tiijj = 3 .* F1./F2;
Tijij = Tiijj./3 + 2./F3 + 1./F4 + (F4.*F5 + F6.*F7 - F8.*F9)./(F2.*F4);
P = Tiijj/3;
Q = (Tijij - P)/5;

% elastic moduli
Ksat=((Phi.*(Kinc-Kmat).*P).*4/3.*Gmat+Kmat.*(Kmat+4/3.*Gmat))./(Kmat+4/3.*Gmat-(Phi.*(Kinc-Kmat).*P));
psi   = (Gmat.*(9*Kmat+8.*Gmat))./(6*(Kmat+2*Gmat));
Gsat= (psi.*(Phi.*(Ginc-Gmat).*Q)+Gmat.*(Gmat+psi))./(Gmat+psi-(Phi.*(Ginc-Gmat).*Q));

% velocities
Vp=sqrt((Ksat+4/3*Gsat)./Rho);
Vs=sqrt(Gsat./Rho);
