function  Vp = WyllieModel(Phi, Vpmat, Vpfl)

% WYLLIE MODEL implements Wyllie's equation 
% INPUT Phi = Porosity
%       Vpmat = P-wave velocity of the solid phase
%       Vpfl = P-wave velocity of the fluid phase
% OUTUPT Vp = P-wave velocity of saturated rock

% Written by Dario Grana (August 2020)

% Wyllie 
Vp = 1./((1-Phi)./Vpmat+Phi./Vpfl);
