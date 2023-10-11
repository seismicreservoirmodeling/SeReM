function [correlation_function] = construct_correlation_function(signal, l, type, angles, taper)
% CORRELATION FUNCTION compute a matrix according to the correlation function, ranges and angles
% Obs: of course it is possible to implement this function using meshgrids
% instead of for loops. However, for large dimension models, meshgrid uses
% too much RAM memory.
% INPUT
%    signal - It is just to define the size of the simulations (FFTMA works better if we define the filter with the same size of the white noise)
%    l - Correlation range, (1,1) or (1, 3)
%    type = function ype ('exp', 'gau', 'sph')
%    angles - used for anisotropic correlation function (1, 3)
%    taper - boolean to apply, or not, the window function. it avoinds artifacts when correlation ranges are similar to the grid dimensions.

I = size(signal,1);
J= size(signal,2);
K = size(signal,3);

% Taper parameters 
desvio = 0.5;
order = 4;

ROT = rotz(angles(3))*roty(angles(2))*rotx(angles(1));
correlation_function = zeros(I,J,K);
for i=1:I
    for j=1:J
        for k=1:K
                                                
            x = (i-round(I/2));
            y = (j-round(J/2));
            z = (k-round(K/2));  
              
            coord = [x;y;z];        
            coord_rot = ROT * coord ;
            coord_rot = coord_rot./l';
            x = coord_rot(1);
            y = coord_rot(2);
            z = coord_rot(3);
            
            h = sqrt( x.^2 + y.^2 + z.^2 );
            
            value = SpatialCovariance1D(h, 1, type);
            
            % taper to avoid FFT artefacts:
            if nargin > 4
                value_window = exp( -(abs((i-round(I/2))/(desvio*I))^order + abs((j-round(J/2))/(desvio*J))^order  + abs((k-round(K/2))/(desvio*K))^order));            
                correlation_function(i,j,k) = value*value_window;
            else
                correlation_function(i,j,k) = value;
            end                        
            
        end
    end
end

correlation_function(round(I/2),round(J/2),round(K/2)) = 1;

