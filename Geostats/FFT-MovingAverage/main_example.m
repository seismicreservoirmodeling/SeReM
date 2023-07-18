

%% Important: the simulations performed by the FFT_MA_3D function are periodic due to the periodic assumption of FFT.
% The usual way to overcome this fact is to generate a simulation larger
% than your model grid and then crop it.
% See Example 2 below as an example. 



%%  EXAMPLE 1
% Variogram types

% Model size:
I = 100;
J = 100;

noise = randn(I,J);

[correlation_function_exp] = construct_correlation_function(10, 30, noise, 'exp', 0);
[ simulation_exp ] = FFT_MA_3D( correlation_function_exp, noise );

[correlation_function_gau] = construct_correlation_function(10, 30, noise, 'gau', 30);
[ simulation_gau ] = FFT_MA_3D( correlation_function_gau, noise );

[correlation_function_sph] = construct_correlation_function(10, 30, noise, 'sph', 60);
[ simulation_sph ] = FFT_MA_3D( correlation_function_sph, noise );


figure
subplot(131)
imagesc(simulation_exp)
subplot(132)
imagesc(simulation_gau)
subplot(133)
imagesc(simulation_sph)

%%  EXAMPLE 2
% Local variable mean and std
% The simulations performed by the FFT_MA_3D function are periodic due to the periodic assumption of FFT.
% The usual way to overcome this fact is to generate a simulation larger than your model grid and then crop it.


% Increasing the grid model to avoid periodicity
I = 2 * I;
J = 2 * J;

noise = randn(I,J);
[correlation_function_sph] = construct_correlation_function(10, 10, noise, 'gau', 0);
[ simulation_sph ] = FFT_MA_3D( correlation_function_sph, noise );

% croping the simulation to avoid periodicity
simulation_sph = simulation_sph(1:I/2,1:J/2);


figure
imagesc(simulation_sph )





