function [Kmat, Gmat, Rhomat, Kfl, Rhofl] = MatrixFluidModel (Kminc, Gminc, Rhominc, Volminc, Kflc, Rhoflc, Sflc, patchy)

% MATRIX FLUID MODEL computes elastic moduli and density of the solid phase
% and fluid phase using Voigt-Reuss averages
% INPUT Kminc = Row vector of bulk moduli of minerals in GPa (ex [36 21])
%       Gminc = Row vector of shear moduli of minerals in GPa (ex [45 7])
%       Rhominc = Row vector of densities of minerals in g/cc (ex [2.6 2.3])
%       Volminc = Matrix of volumes. Each column is a mineral volume log
%                 (ex [vquartz 1-vclay])
%       Kflc = Row vector of bulk moduli of fluid components in GPa 
%               (ex [2.25 0.8 0.1])
%       Rhoflc = Row vector of densities of fluid components in g/cc 
%               (ex [1.03 0.7 0.02])
%       Volminc = Matrix of saturations. Each column is a saturation log
%                 (ex [sw so sg 1-vclay])
%       patchy = binary variable: 1=Patchy; 0=Homegeneous
% OUTUPT Kmat = bulk modulus of matrix phase
%        Gmat = shear modulus of matrix phase
%        Rhomat = density of matrix phase
%        Kfl = bulk modulus of fluid phase
%        Rhofl = density of fluid phase

% Written by Dario Grana (August 2020)

% number of samples
n=size(Volminc,1);
% initialization variables
KmatV=zeros(n,1);
KmatR=KmatV; Kmat=KmatV;
GmatV=KmatV; GmatR=KmatV; Gmat=KmatV;
Rhomat=KmatV; Kfl=KmatV; Rhofl=KmatV;
for i=1:size(Volminc,1)
    % Voigt average (bulk)
    KmatV(i)=sum((Volminc(i,:).*Kminc)/sum(Volminc(i,:)));
    % Reuss average (bulk)
    KmatR(i)=1./sum((Volminc(i,:)./Kminc)/sum(Volminc(i,:)));
    % Voigt-Reuss-Hill average (bulk)
    Kmat(i)=0.5*(KmatV(i)+KmatR(i));
    % Voigt average (shear)
    GmatV(i)=sum((Volminc(i,:).*Gminc)/sum(Volminc(i,:)));
    % Reuss average (shear)
    GmatR(i)=1./sum((Volminc(i,:)./Gminc)/sum(Volminc(i,:)));
    % Voigt-Reuss-Hill average (shear)
    Gmat(i)=0.5*(GmatV(i)+GmatR(i));
    % linear average for matrix density
    Rhomat(i)=sum((Volminc(i,:).*Rhominc)/sum(Volminc(i,:)));
    if patchy==0
        % Reuss average for fluid
        Kfl(i)=1./sum(Sflc(i,:)./Kflc);
    else
        % Voigt average for fluid
        Kfl(i)=sum(Sflc(i,:).*Kflc);
    end
    % linear average for fluid density
    Rhofl(i)=sum(Sflc(i,:).*Rhoflc);
end

