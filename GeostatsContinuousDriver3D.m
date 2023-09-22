
addpath(genpath('../SeReM/'))
close all

%% Example 1
% available data (4 measurements)
%dcoords = [ 5 18; 15 13; 11 4; 1 9 ];
dcoords = [ 5 18 0; 15 13 0; 11 4 0; 1 9 0];

dvalues = [ 3.1 3.9 4.1 3.2 ]';
% coordinates of the location to be estimated
xcoords = [ 11 9 0];

% parameters random variable
xmean = 3.5;
xvar = 0.1;
l = [90 2 90];
%l = [2 2 2];
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

%% Example 
 load Data/ElevationData.mat
 % available data (15 measurements)

Lx = 18;
Ly = 20;
Lz = 22;
dcoords = [ Ly/4   Lx/2   Lz/2 ;
           3*Ly/4  Lx/2  Lz/2 ];        
dcoords = round(dcoords);
dz = [ -1 ; 1 ];

%[X,Y,T] = meshgrid( 1:8, 1:10, 1:12 );
[Y,X,T] = meshgrid( 1:Lx, 1:Ly, 1:Lz );
xcoords = [ X(:) Y(:) T(:) ];
n = size(xcoords,1);

% parameters random variable
zmean = 0;
zvar = 1;
lx = 15;
ly = 15;
lz = 15;
l = [ly lx lz];
type = 'gau';


% kriging                    
[xok_grid, v_grid] = Kriging_options(xcoords, dcoords, dz, zvar, l, type);
xok_grid = reshape(xok_grid, size(X));
v_grid = reshape(v_grid, size(X));
v_grid(v_grid<0) = 0;

% % Sequential Gaussian Simulation
% krig = 1;
% sgsim = SeqGaussianSimulation(xcoords, dcoords, dz, zmean, zvar, l, type, krig);
% sgsim = reshape(sgsim,size(X));
   

% plot results
figure(3)
subplot(3,3,1)
imagesc(squeeze(xok_grid(:,:,end/2)))
caxis([-2.5 2.5])
xlabel('X')
ylabel('Y')
subplot(3,3,2)
imagesc(squeeze(xok_grid(:,end/2,:)))
caxis([-2.5 2.5])
xlabel('Z')
ylabel('Y')
subplot(3,3,3)
imagesc(squeeze(xok_grid(end/2-1,:,:)))
caxis([-2.5 2.5])
xlabel('Z')
ylabel('X')
subplot(3,3,4)
imagesc(squeeze(v_grid(:,:,end/2)))
caxis([0 1])
xlabel('X')
ylabel('Y')
subplot(3,3,5)
imagesc(squeeze(v_grid(:,end/2,:)))
caxis([0 1])
xlabel('Z')
ylabel('Y')
subplot(3,3,6)
imagesc(squeeze(v_grid(end/2,:,:)))
caxis([0 1])
xlabel('Z')
ylabel('X')
% subplot(3,3,7)
% imagesc(squeeze(sgsim(:,:,end/2)))
% caxis([-2.5 2.5])
% xlabel('X')
% ylabel('Y')
% subplot(3,3,8)
% imagesc(squeeze(sgsim(:,end/2,:)))
% caxis([-2.5 2.5])
% xlabel('Z')
% ylabel('Y')
% subplot(3,3,9)
% imagesc(squeeze(sgsim(end/2,:,:)))
% caxis([-2.5 2.5])
% xlabel('Z')
% ylabel('X')


%% Example 2 (elevation Yellowstone)
% load Data/ElevationData.mat
% % available data (15 measurements)
% dt = zeros(size(dy));
% %dcoords = [ dx dy dt]; % for the set of 100 measurements simply comment line 56
% dcoords = [ dy dx dt]; % for the set of 100 measurements simply comment line 56
% nd = size(dcoords,1);
% % grid of coordinates of the location to be estimated
% %[Y,X,T] = meshgrid( X(1,:), Y(:,1), 0:2 );
% [Y,X,T] = meshgrid(X(1,:), Y(:,1), 0:2 );
% xcoords = [ X(:) Y(:) T(:) ];
% n = size(xcoords,1);
% 
% % parameters random variable
% zmean = 2476;
% zvar = 8721;
% lx = 12.5;
% ly = 6;
% lz = 12.5;
% l = [ly lx lz];
% type = 'exp';
% 
% % plot
% figure(4)
% scatter(dx, dy, 50, dz, 'filled');
% grid on; box on; xlabel('X'); ylabel('Y'); colorbar; caxis([2000 2800]); 
% 
% % kriging
%  [xok_grid, v_grid] = Kriging_options(xcoords, dcoords, dz, zvar, l, type);
%  xok_grid = reshape(xok_grid,size(X));
%  v_grid = reshape(v_grid,size(X));
% 
% % Sequential Gaussian Simulation
% krig = 1;
% sgsim = SeqGaussianSimulation(xcoords, dcoords, dz, zmean, zvar, l, type, krig);
% sgsim = reshape(sgsim,size(X));
% 
% 
% % plot results
% figure(4)
% subplot(221)
% pcolor(unique(Y),unique(X),xok_grid(:,:,1));
% xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
% colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
% title('Simple Kriging');
% subplot(222)
% imagesc(unique(Y),unique(X),xok_grid(:,:,3));
% xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
% colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
% title('Ordinary Kriging');
% subplot(223)
% imagesc(unique(Y),unique(X),sgsim(:,:,1));
% xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
% colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
% title('SGS Realization 1');
% subplot(224)
% imagesc(unique(Y),unique(X),sgsim(:,:,3));
% xlabel('X'); ylabel('Y'); shading interp; set(gca, 'YDir', 'reverse')
% colorbar; caxis([2000 2800]); hbc=colorbar; title(hbc, 'Elevation');
% title('SGS Realization 2');








