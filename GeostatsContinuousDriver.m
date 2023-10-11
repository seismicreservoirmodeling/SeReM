%% Geostatistics Continuous Driver %%
% In this script we illustrate kriging and sequential simulation with 
% two examples: 
% Example 1: example with 4 density measurements
% Example 2: example with 15 elevation measurements from Yellowstone

addpath(genpath('../SeReM/'))

%% Example 1
% available data (4 measurements)
dcoords = [ 5 18; 15 13; 11 4; 1 9 ];
dvalues = [ 3.1 3.9 4.1 3.2 ]';
% coordinates of the location to be estimated
xcoords = [ 10 10 ];

% parameters random variable
xmean = 3.5;
xvar = 0.1;
l = 9;
type = 'exp';

% plot
figure(1)
scatter(dcoords(:,1), dcoords(:,2), 100, dvalues, 'filled');
hold on
plot(xcoords(:,1), xcoords(:,2), 'ks');
grid on; box on; xlabel('X'); ylabel('Y'); colorbar; caxis([3 4.2]); 

% simple kiging
[xsk, ~] = SimpleKriging(xcoords, dcoords, dvalues, xmean, xvar, l, type, [0 0 0]);

% ordinary kiging
[xok, ~] = OrdinaryKriging(xcoords, dcoords, dvalues, xvar, l, type, [0 0 0]);

% Gaussian simulation
krig = 0;
nsim = 100;
gsim = zeros(nsim,1);
for i=1:nsim
    gsim(i) = GaussianSimulation(xcoords, dcoords, dvalues, xmean, xvar, l, type, krig, [0 0 0]);
end

% plot results
figure(2)
histogram(gsim);
hold on
plot(xsk, 0, '*r')
plot(xok, 0, 'sb')
plot(mean(gsim), 0, 'og')
grid on; box on; xlabel('Property'); ylabel('Frequency'); 
legend('Gauss Sims.', 'Simple Krig.', 'Ord. Krig.', 'mean Gauss Sims.');


%% Example 2 (elevation Yellowstone)
load Data/ElevationData.mat
% available data (15 measurements)
dcoords = [ dx dy ]; % for the set of 100 measurements simply comment line 56
nd = size(dcoords,1);
% grid of coordinates of the location to be estimated
xcoords = [ X(:) Y(:) ];
n = size(xcoords,1);

% parameters random variable
zmean = 2476;
zvar = 8721;
l = 12.5;
type = 'exp';

% plot
figure(3)
scatter(dcoords(:,1), dcoords(:,2), 50, dz, 'filled');
grid on; box on; xlabel('X'); ylabel('Y'); colorbar; caxis([2000 2800]); 

% kriging
xsk = zeros(nd,1);
xok = zeros(nd,1);
for i=1:n
    % simple kiging
    [xsk(i), ~] = SimpleKriging(xcoords(i,:), dcoords, dz, zmean, zvar, l, type);
    % ordinary kiging
    [xok(i), ~] = OrdinaryKriging(xcoords(i,:), dcoords, dz, zvar, l, type);
end
xsk = reshape(xsk,size(X));
xok = reshape(xok,size(X));

% Sequential Gaussian Simulation
krig = 1;
nsim = 3;
sgsim = zeros(size(X,1),size(X,2),nsim);
for i=1:nsim
    sim = SeqGaussianSimulation(xcoords, dcoords, dz, zmean, zvar, l, type, krig, [0 0 0]);
    sgsim(:,:,i) = reshape(sim,size(X,1),size(X,2));
end

% plot results
figure(4)
subplot(221)
pcolor(unique(X),unique(Y),xsk);
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
title('Simple Kriging');
subplot(222)
imagesc(unique(X),unique(Y),xok);
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
title('Ordinary Kriging');
subplot(223)
imagesc(unique(X),unique(Y),sgsim(:,:,1));
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
title('SGS Realization 1');
subplot(224)
imagesc(unique(X),unique(Y),sgsim(:,:,2));
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
title('SGS Realization 2');


%% Example 3 (elevation Yellowstone) using grid Kriging_options function
% Kriging_options is being implemented to account for additional options
% such as max number of conditioning points, searching neighborhood and type
% of kriging, including the parfor loop over all the points 

% resampling random positions to test  max number of conditioning points
n_dcoords = 200;
sampling_indices = randperm(numel(X(:)), n_dcoords );
dcoords = [ X(sampling_indices)' Y(sampling_indices)' ];
dz = Z(sampling_indices)';
xcoords = [ X(:) Y(:) ];
n = size(xcoords,1);
nd = size(dcoords,1);

% parameters random variable
zmean = 2476;
zvar = 8721;
l = 12.5;
type = 'exp';

% plot
figure(5)
subplot(131)
scatter(dcoords(:,1), dcoords(:,2), 50, dz, 'filled');
grid on; box on; xlabel('X'); ylabel('Y'); colorbar; caxis([2000 2800]); 
set(gca, 'YDir', 'reverse')
axis([min(dcoords(:,1)) max(dcoords(:,1)) min(dcoords(:,2)) max(dcoords(:,2)) ])

[xok_grid, v_grid] = Kriging_options(xcoords, dcoords, dz, zvar, l, type, krig, [0 0 0]);

xok_grid = reshape(xok_grid, size(X));
v_grid = reshape(v_grid, size(X));
v_grid(v_grid<0) = 0;

subplot(132)
pcolor(unique(X),unique(Y),xok_grid);
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
title('Result of Kriging-options')

subplot(133)
pcolor(unique(X),unique(Y),v_grid);
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar;  hbc=colorbar; title(hbc, 'Elevation');
title('Result of Kriging-options')


%% Example 4 FFT-Moving average simulation example
%% Important: the simulations performed by the FFT_MA_3D function are periodic due to the periodic assumption of FFT.
% The usual way to overcome this fact is to generate a simulation larger than your model grid and then crop it.

I = 2 * size(X,1);
J = 2 * size(X,2);

noise = randn(I,J);

[correlation_function] = construct_correlation_function(noise, 15, type, [0 0 0]);
[ simulation ] = FFT_MA_3D( correlation_function, noise );

% croping the simulation to avoid periodicity
simulation = simulation(1:I/2,1:J/2);

% Based on PFS, apply kriging mean and variance for conditional simulations
simulation = xok_grid + sqrt(v_grid) .* simulation;


% plot
figure(6)
pcolor(unique(X),unique(Y),simulation);
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');




%% Example 5 DMS - Direct Multivariate Simulation example

load('Data/HardData_ReferenceModel_size100_range20.mat');

% Use 2D reference simulations as reference variables:
I = size(reference_models,2);
J = size(reference_models,3);
reference_variables = [reshape(reference_models(1,:,:),I*J,1) reshape(reference_models(2,:,:),I*J,1) reshape(reference_models(3,:,:),I*J,1) reshape(reference_models(4,:,:),I*J,1) reshape(reference_models(5,:,:),I*J,1) reshape(reference_models(6,:,:),I*J,1) ] ;

% DKE: OPTIONAL
[reference_variables] = extend_dateset_KDE(reference_variables,2,0.05);

% Run DMS 
cell_size = 0.05; 
l = 20;
n_simulations = 1;

% Number of conditional points
n_cond_points = 100;
cond_value_ = cond_value(1:n_cond_points ,:);
cond_pos_ = cond_pos(1:n_cond_points ,:);

% condicional DMS
[simulations_all_dms] = DMS(I,J, l, type, [0 0 0], cell_size, reference_variables, cond_pos_, cond_value_, n_simulations);
simulation_dms = simulations_all_dms{1};


% plot
% simulations
generate_2D(reference_models,cond_pos_)
subplot(2,3,1)
caxis([-3 3])
subplot(2,3,2)
caxis([-15 2])
subplot(2,3,3)
caxis([-10 12])
subplot(2,3,4)
caxis([-3 3])
subplot(2,3,5)
caxis([-5 5])
subplot(2,3,5)
caxis([-13 10])

generate_2D(simulation_dms,cond_pos_)
subplot(2,3,1)
caxis([-3 3])
subplot(2,3,2)
caxis([-15 2])
subplot(2,3,3)
caxis([-10 12])
subplot(2,3,4)
caxis([-3 3])
subplot(2,3,5)
caxis([-5 5])
subplot(2,3,5)
caxis([-13 10])

% Histograms
generate_histograms(reshape(reference_models,6,I*J)')
generate_histograms(reshape(simulation_dms,6,I*J)')
















