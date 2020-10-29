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
[xsk, ~] = SimpleKriging(xcoords, dcoords, dvalues, xmean, xvar, l, type);

% ordinary kiging
[xok, ~] = OrdinaryKriging(xcoords, dcoords, dvalues, xvar, l, type);

% Gaussian simulation
krig = 0;
nsim = 100;
gsim = zeros(nsim,1);
for i=1:nsim
    gsim(i) = GaussianSimulation(xcoords, dcoords, dvalues, xmean, xvar, l, type, krig);
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
load Data/data6.mat
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
    sim = SeqGaussianSimulation(xcoords, dcoords, dz, zmean, zvar, l, type, krig);
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

