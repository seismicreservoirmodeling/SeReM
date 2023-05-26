function [Vp, Vs, Rho] = RPM(Phi ,v_clay, sw)

criticalporo=0.4;
coordnumber=9;
pressure=0.060;


Kminc = [37 21];
Gminc = [45 7];
Rhominc = [2.65 2.58];
Volminc = [ 1-v_clay  v_clay ] ;

Kflc = [3.14 0.53];
Rhoflc = [1.06 0.52];
Sflc = [sw 1-sw];

[Kmat, Gmat, Rhomat, Kfl, Rhofl] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, 0);

Rho = DensityModel(Phi, Rhomat, Rhofl);

[Vp, Vs] = SoftsandModel(Phi, Rho, Kmat, Gmat, Kfl, criticalporo, coordnumber, pressure);