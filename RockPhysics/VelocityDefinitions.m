function  [Vp, Vs] = VelocityDefinitions(K, G, Rho)

% VELOCITY DEFINITIONS implements the definitions of P- and S-wave velocity
% INPUT K = Bulk modulus
%       G = Shear modulus
%       Rho = Density
% OUTUPT Vp = P-wave velocity
%        Vs = S-wave velocity

% Written by Dario Grana (August 2020)

% definitions
Vp=sqrt((K+4/3*G)/Rho);
Vs=sqrt(G/Rho);
