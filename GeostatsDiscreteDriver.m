%% Geostatistics Discrete Driver %%
% In this script we illustrate kriging and sequential simulation with 
% two examples: 
% Example 1: example with 4 density measurements
% Example 2: example with 15 elevation measurements from Yellowstone

addpath(genpath('../SeReM/'))

%% Example 1
% available data (4 measurements)
dcoords = [ 5 18; 15 13; 11 4; 1 9 ];
fvalues = [ 1 2 2 1 ]';
% coordinates of the location to be estimated
xcoords = [ 10 10 ];

% parameters random variable
nf = 2;
pprior = [ 0.5 0.5];
l = 9;
type = 'exp';

% plot
figure(1)
scatter(dcoords(:,1), dcoords(:,2), 100, fvalues, 'filled');
hold on
plot(xcoords(:,1), xcoords(:,2), 'ks');
grid on; box on; xlabel('X'); ylabel('Y'); colorbar; 

% indicator kriging
[ikp, ikmap] = IndicatorKriging(xcoords, dcoords, fvalues, nf, pprior, l, type);

% simulation
nsim = 1000;
isim = zeros(nsim,1);
for i=1:nsim
     isim(i) = RandDisc(ikp);
end

% plot results
figure(2)
histogram(isim, 'BinWidth', 0.1);
hold on
plot(ikmap, 0, '*r')
grid on; box on; xlabel('Discrete property'); ylabel('Frequency'); 
legend('Ind. Krig.', 'Ind. Krig. MAP');

%% Example 2 (elevation Yellowstone)
load Data/ElevationData.mat
load Data/data6.mat
% available data (15 measurements)
dcoords = [ dx dy ]; % for the set of 100 measurements simply comment line 56
nd = size(dcoords,1);

% discrete property definition
zmean = 2476;
df = ones(size(dz));
df(dz>zmean) = 2;

% grid of coordinates of the location to be estimated
xcoords = [ X(:) Y(:) ];
n = size(xcoords,1);

% parameters random variable
pprior = [ 0.5 0.5 ];
l = 12.5;
type = 'exp';

% plot
figure(3)
scatter(dcoords(:,1), dcoords(:,2), 50, df, 'filled');
grid on; box on; xlabel('X'); ylabel('Y'); colorbar; 

% kriging
ikp = zeros(n,nf);
for i=1:n
    [ikp(i,:), ikmap(i)] = IndicatorKriging(xcoords(i,:), dcoords, df, nf, pprior, l, type);
end
ikp = reshape(ikp,size(X,1),size(X,2),nf);
ikmap = reshape(ikmap,size(X,1),size(X,2));

% Sequential Indicator Simulation
nsim = 3;
sisim = zeros(size(X,1),size(X,2),nsim);
for i=1:nsim
    sim = SeqIndicatorSimulation(xcoords, dcoords, df, nf, pprior, l, type);
    sisim(:,:,i) = reshape(sim,size(X,1),size(X,2));
end

% plot results
figure(4)
subplot(221)
pcolor(unique(X),unique(Y),ikp(:,:,1));
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; caxis([0 1]); 
title('Indicator Kriging Probability');
subplot(222)
imagesc(unique(X),unique(Y),ikmap);
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; 
title('Indicator Kriging MAP');
subplot(223)
imagesc(unique(X),unique(Y),sisim(:,:,1));
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; title('SIS Realization 1');
subplot(224)
imagesc(unique(X),unique(Y),sisim(:,:,2));
xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
colorbar; title('SIS Realization 2');


%% Markov chain simulation
% initial parameters
nsim = 3;
ns = 100;
% vertical axis
z = (ns:-1:1)';

% Transition matrix T1 (equal propotions, equal transitions)
T1 = [ 0.5 0.5
       0.5 0.5];
% Transition matrix T2 (equal propotions, asymmetrix transitions)
T2 = [ 0.9 0.1
       0.1 0.9];
% Transition matrix T3 (different propotions, asymmetrix transitions)
T3 = [ 0.1 0.9
       0.1 0.9];
   
% simulation
fsim1 = MarkovChainSimulation(T1, ns, nsim);
fsim2 = MarkovChainSimulation(T2, ns, nsim);
fsim3 = MarkovChainSimulation(T3, ns, nsim);

% plot realzations
figure(5)
subplot(131)
imagesc((1:nsim),z,fsim1);
xlabel('Facies realizations'); ylabel('Relative depth (m)'); box on;
title('Transition matrix T_1')
subplot(132)
imagesc((1:nsim),z,fsim2);
xlabel('Facies realizations'); set(gca, 'YTickLabel', []); box on;
title('Transition matrix T_2')
subplot(133)
imagesc((1:nsim),z,fsim3);
xlabel('Facies realizations'); set(gca, 'YTickLabel', []); box on;
title('Transition matrix T_3')
