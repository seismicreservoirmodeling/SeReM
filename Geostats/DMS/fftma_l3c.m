function [ noises ] = fftma_l3c(m,n, range_h, range_v,err, root_noises)
%FFTMA_L3C Summary of this function goes here
%   Detailed explanation goes here
%FFTMA
range_vertical = range_v;
range_horizontal = range_h;
noises = root_noises;
noises = reshape(noises,m,n);
[correlation_function] = construct_correlation_function(range_vertical, range_horizontal, noises, 'sph', 0);        

noises = FFT_MA_3D(correlation_function,noises);

noises = reshape(noises, m , n);

%noises = noises-mean(noises(:));
%noises = noises./std(noises(:));
noises = make_it_gaussian(noises(:));

noises = reshape(noises,m,n);

end

