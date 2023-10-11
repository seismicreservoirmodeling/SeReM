function [correlation_function] = construct_correlation_function(signal, l, type, angles)
% Obs: of course it is possible to implement this function using meshgrids
% instead of for loops. However, for large dimension models, meshgrid uses
% too much RAM memory.
% Lv - Vertical correlation range
% Lh - Horizontal correlation range
% signal - It is just to define the size of the simulations (FFTMA works better if we define the filter with the same size of the white noise)
% type - 1 for Gaussian, 1 for exponential and 3 for spherical
% angle - Used for anisotropic spherical models

I = size(signal,1);
J= size(signal,2);
K = size(signal,3);

% Taper parameters (it avoinds artifacts when the variogram range are similar to the model size):
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
            %value_window = exp( -(abs((i-round(I/2))/(desvio*I))^order + abs((j-round(J/2))/(desvio*J))^order  + abs((k-round(K/2))/(desvio*K))^order));            
            %correlation_function(i,j,k) = value*value_window;
            
            correlation_function(i,j,k) = value;
            
        end
    end
end

correlation_function(round(I/2),round(J/2),round(K/2)) = 1;

