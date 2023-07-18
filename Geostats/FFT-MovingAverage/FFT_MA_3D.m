function [ simulation ] = FFT_MA_3D( correlation_function, noise )


%% Treating the signal to be power of 2
%  S = size(noise);
%  S = log2(S);
%  S = round(S+0.5);
%  S = 2.^S;
%  
%  c2 = ifftn(sqrt(abs(fftn(correlation_function, S ))) .* fftn(noise,S),S);
%  
%  simulation = real(c2);
%  simulation = simulation(1:size(noise,1), 1:size(noise,2));

%% Atraight forward
c2 = ifftn(sqrt(abs(fftn(correlation_function,size(noise)))) .* fftn(noise,size(noise)));

simulation = real(c2);


end