function  Vp = RaymerModel(Phi, Vpmat, Vpfl)

% RAYMER MODEL implements Raymer's equation 
% INPUT Phi = Porosity
%       Vpmat = P-wave velocity of the solid phase
%       Vpfl = P-wave velocity of the fluid phase
% OUTUPT Vp = P-wave velocity of saturated rock

% Written by Dario Grana (August 2020)

% Raymer  
Vp = (1-Phi).^2*Vpmat+Phi.*Vpfl;


