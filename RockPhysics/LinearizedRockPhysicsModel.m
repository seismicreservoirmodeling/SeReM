function [Vp, Vs, Rho] = LinearizedRockPhysicsModel(Phi, Clay, Sw, R)

% LINEARIZED ROCK PHYSICS MODEL implements a linear rock physics model
% based on a multilinear regression 
% INPUT Phi = Porosity
%       Clay = Clay volume
%       Sw = Shear modulus of dry rock
%       R = regression coefficients matrix (estimated with regress.m)
% OUTUPT Vp = P-wave velocity
%        Vs = S-wave velocity
%        Rho = Density

% Written by Dario Grana (August 2020)

% multilinear regression
Vp = R(1,1)*Phi+R(1,2)*Clay+R(1,3)*Sw+R(1,4);
Vs = R(2,1)*Phi+R(2,2)*Clay+R(2,3)*Sw+R(2,4);
Rho = R(3,1)*Phi+R(3,2)*Clay+R(3,3)*Sw+R(3,4);