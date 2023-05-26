function [ simulation ] = FFT_MA_3D( correlation_function, noise )

	X = correlation_function;
	Y = noise;

    c2 = ifftn(sqrt(abs(fftn(X,size(Y)))) .* fftn(Y,size(Y)));

    simulation = real(c2);
    

end