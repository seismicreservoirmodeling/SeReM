function [corrfun] = CorrelationFunction3D(Lv,Lh,nxl,nil,nm)

% CorrelationFunction3D computes the 3D spatial correlation function 
% using FFT-MA simulations
% INPUT vertcorr = vertical correlation parameter
%       horcorr = horizontal correlation parameter
%       nxl = number of crossline
%       nil - number of inline
%       nm = number of samples
% OUTPUT corrfun = Spatial Correlation model

% Written by Leandro de Figueiredo (March 2018)
% Modified by Dario Grana (June 2023)

ordem = 3;
desvio = 0.25;
corrfun = zeros(nxl,nil,nm);
for i=1:nxl
    for j=1:nil
        for k=1:nm
            r = sqrt( ((i-round(nxl/2))/(3*Lh))^2 + ((j-round(nil/2))/(3*Lh))^2 + ((k-round(nm/2))/(3*Lv))^2);
            if r<1
                value = 1 - 1.5 * r + 0.5 * r^3;
            else
                value=0;
            end
            winval = exp(-(abs((i-round(nxl/2))/(desvio*nxl))^ordem + abs((j-round(nil/2))/(desvio*nil))^ordem  + abs((k-round(nm/2))/(desvio*nm))^ordem));
            corrfun(i,j,k) = value*winval;
        end
    end
end

