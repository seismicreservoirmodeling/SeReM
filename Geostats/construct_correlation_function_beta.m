function [correlation_function] = construct_correlation_function_beta(Lv, Lh, signal, type)
% Lv - Vertical correlation range
% Lh - Horizontal correlation range
% signal - It is just to define the size of the simulations (FFTMA works better if we define the filter with the same size of the white noise)

ordem = 4;
I = size(signal,1);
J= size(signal,2);
K = size(signal,3);
desvio = 1.0;

correlation_function = zeros(I,J,K);
for i=1:I
    for j=1:J
        for k=1:K
             
            if type==1
                %value = exp( -((((i-round(I/2))^2)/Lv^2) + (((j-round(J/2))^2)/Lh^2) + (((k-round(K/2))^2)/Lh^2) )); % )); % 
                value = exp( -(sqrt(((i-round(I/2))^2)/Lv^2) + sqrt(((j-round(J/2))^2)/Lh^2) + sqrt(((k-round(K/2))^2)/Lh^2) ));
            
            else
                r = sqrt( ((i-round(I/2))/(3*Lv))^2 + ((j-round(J/2))/(3*Lh))^2 + ((k-round(K/2))/(3*2*Lh))^2);
                if r<1
                    value = 1 - 1.5 * r + 0.5 * r^3;
                else
                    value=0;
                end
            end
            
            value_window = exp( -(abs((i-round(I/2))/(desvio*I))^ordem + abs((j-round(J/2))/(desvio*J))^ordem  + abs((k-round(K/2))/(desvio*K))^ordem));
            correlation_function(i,j,k) = value*value_window;
            %correlation_function(i,j,k) = value;
            
        end
    end
end

