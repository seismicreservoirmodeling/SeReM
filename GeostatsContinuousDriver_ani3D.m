
%% OBSERVATION: The angles of the anisotropic variogram are still hard-coded in Ordinary and Simples Kriging functions

addpath(genpath('../SeReM/'))

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
angles = [0 0 45];
%l = [2 2 2];
type = 'exp';

% plot
figure(1)
scatter(dcoords(:,1), dcoords(:,2), 100, dvalues, 'filled');
hold on
plot(xcoords(:,1), xcoords(:,2), 'ks');
grid on; box on; xlabel('X'); ylabel('Y'); colorbar; caxis([3 4.2]); 

% simple kiging
[xsk, ~] = SimpleKriging(xcoords, dcoords, dvalues, xmean, xvar, l, type, angles);

% ordinary kiging
[xok, ~] = OrdinaryKriging(xcoords, dcoords, dvalues, xvar, l, type, angles);

% Gaussian simulation
krig = 1;
nsim = 100;
gsim = zeros(nsim,1);
for i=1:nsim
    gsim(i) = GaussianSimulation(xcoords, dcoords, dvalues, xmean, xvar, l, type, krig, angles);
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

% grid size:
Li = 20;
Lj = 20;
Lk = 20;
% Matlab meshgrid invert the dimension of X and Y when working with 3
% dimension, because of that we need to invert them in the following lines:
[J,I,K] = meshgrid( 1:Li, 1:Lj, 1:Lk );
xcoords = [ I(:) J(:) K(:) ];

% Hard data:
i_data = [Li/2 ; Li/2];
j_data = [Lj/4 ; 3*Lj/4];
k_data = [Lk/2 ; Lk/2];
% We need to invert x y here as well
dcoords = [j_data i_data k_data];
dcoords = round(dcoords);
dz = [ -1 ; 1 ];

% Random variable / Variogram settings
zmean = 0;
zvar = 1;
Li = 5;
lj = 15;
lk = 5;
% We need to invert x y here as well
l = [lj Li lk];

angi = 45;
angj = 0;
angk = 45;
% We need to invert x y here as well
angles = [angj angi angk];

type = 'gau';
krig = 1;

% kriging                    
[krig_mean, krig_var] = Kriging_options(xcoords, dcoords, dz, zvar, l, type, krig, angles);
krig_mean = reshape(krig_mean, size(I));
krig_mean = permute(krig_mean, [2, 1, 3]);
krig_var = reshape(krig_var, size(I));
krig_var(krig_var<0) = 0;
% We need to invert x y here as well
krig_var = permute(krig_var, [2, 1, 3]);

% I do not have any idea why, but the SGS is working only with Li=lj=lk,
% otherwise it gives several outliers.
% % Sequential Gaussian Simulation
krig = 1;
sgsim = SeqGaussianSimulation(xcoords, dcoords, dz, zmean, zvar, l, type, krig, angles);
sgsim = reshape(sgsim,size(I));
% We need to invert x y here as well
sgsim = permute(sgsim , [2, 1, 3]);
   

% plot results
figure(3)
subplot(3,3,1)
imagesc(squeeze(krig_mean(:,:,end/2)))
caxis([-2 2])
xlabel('Y')
ylabel('X')
subplot(3,3,2)
imagesc(squeeze(krig_mean(:,Lj/4,:)))
caxis([-2 2])
xlabel('Z')
ylabel('X')
subplot(3,3,3)
imagesc(squeeze(krig_mean(end/2,:,:)))
caxis([-2 2])
xlabel('Z')
ylabel('Y')
subplot(3,3,4)
imagesc(squeeze(krig_var(:,:,end/2)))
caxis([0 1])
xlabel('Y')
ylabel('X')
subplot(3,3,5)
imagesc(squeeze(krig_var(:,Lj/4,:)))
caxis([0 1])
xlabel('Z')
ylabel('X')
subplot(3,3,6)
imagesc(squeeze(krig_var(end/2,:,:)))
caxis([0 1])
xlabel('Z')
ylabel('Y')
subplot(3,3,7)
imagesc(squeeze(sgsim(:,:,end/2)))
caxis([-2.5 2.5])
xlabel('Y')
ylabel('X')
subplot(3,3,8)
imagesc(squeeze(sgsim(:,Lj/4,:)))
caxis([-2.5 2.5])
xlabel('Z')
ylabel('X')
subplot(3,3,9)
imagesc(squeeze(sgsim(end/2,:,:)))
caxis([-2.5 2.5])
xlabel('Z')
ylabel('Y')








