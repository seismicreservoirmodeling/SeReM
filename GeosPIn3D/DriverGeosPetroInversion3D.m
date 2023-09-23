clear all

%% Geostatistical Petrophysical inversion Driver %%
% In this script we apply the Ensenmble petrophysical inversion method
% (Liu and Grana, 2018) to predict the petrophysical properties 
% from seismic data.

%% Available data and parameters
% Load data (seismic data and time)
addpath(genpath('../SeReM/'))
addpath(genpath('../GeoSPin3D/'))

% real dataset
load Data3D/SeismicData3D.mat
interx = 39; intery = 34;
% % synthetic dataset
% load Data3D/SeismicData3D_test.mat
% interx = 2; intery = 2;


%% Initial parameters
% number of samples (petrophysical properties)
nm = size(near,3)+1;
% number of lines
nxl = size(near,1);
nil = size(near,2);
% number of variables
nv = 3;
% reflection angles 
theta = [15, 30, 45];
ntheta = length(theta);
% time vector
TimeSeis = squeeze(Z(1,1,:));
% time sampling
dt = TimeSeis(2)-TimeSeis(1);
% error variance
varerr = 10^-3;
sigmaerr = varerr*eye(ntheta*(nm-1));

%% Wavelet
% wavelet 
freq = 45;
ntw = 64;
[wavelet, tw] = RickerWavelet(freq, dt, ntw);

%% Plot seismic data
figure(1)
subplot(131)
h = slice(X,Y,Z,near,interx,intery,[]);
set(h,'edgecolor','none')
zlim([TimeSeis(1) TimeSeis(end)]); xlim([1 nil]); ylim([1 nxl]);
axis tight; grid on; box on; set(gca, 'ZDir', 'reverse')
xlabel('Inline'); ylabel('Crossline'); zlabel('Time (s)'); title('Near')
subplot(132)
h = slice(X,Y,Z,mid,interx,intery,[]);
set(h,'edgecolor','none')
zlim([TimeSeis(1) TimeSeis(end)]); xlim([1 nil]); ylim([1 nxl]);
axis tight; grid on; box on; set(gca, 'ZDir', 'reverse')
xlabel('Inline'); ylabel('Crossline'); zlabel('Time (s)'); title('Mid')
subplot(133)
h = slice(X,Y,Z,far,interx,intery,[]);
set(h,'edgecolor','none')
zlim([TimeSeis(1) TimeSeis(end)]); xlim([1 nil]); ylim([1 nxl]);
axis tight; grid on; box on; set(gca, 'ZDir', 'reverse')
xlabel('Inline'); ylabel('Crossline'); zlabel('Time (s)'); title('Far') 

%% Prior model (filtered well logs)
% nfilt = 3;
% cutofffr = 0.04;
% [b, a] = butter(nfilt, cutofffr);
% Vpprior = filtfilt(b, a, Vp);
% Vsprior = filtfilt(b, a, Vs);
% Rhoprior = filtfilt(b, a, Rho);
phipriormean = 0.2;
claypriormean = 0.23;
swpriormean = 0.6;
phiprior = phipriormean*ones(nxl,nil,nm);
clayprior = claypriormean*ones(nxl,nil,nm);
swprior = swpriormean*ones(nxl,nil,nm);

%% Spatial correlation matrix
corrpetro = [   1   -0.6   -0.5
             -0.6    1    0.2
             -0.5    0.2    1];
stdpetro = [0.035 0.055 0.09];

%% Rock physics parameters
% training dataset
load RockPhysicsTrain.mat
% Xreg = [ones(size(VpTrain)) dtrain];
% breg = regress(PhiTrain,Xreg); 
% PhiTrain = breg(1) + breg(2)*VpTrain + breg(3)*VsTrain + breg(4)*RhoTrain;
% breg = regress(ClayTrain,Xreg); 
% ClayTrain = breg(1) + breg(2)*VpTrain + breg(3)*VsTrain + breg(4)*RhoTrain;
% breg = regress(SwTrain,Xreg); 
% SwTrain = breg(1) + breg(2)*VpTrain + breg(3)*VsTrain + breg(4)*RhoTrain;
petrotrain = [PhiTrain ClayTrain SwTrain];
elastrain = [VpTrain VsTrain RhoTrain];


%% Seismic inversion
% sigmaprior, elastrain, petrotrain, faciestrain, vpgrid, vsgrid, rhogrid, phigrid, claygrid, swgrid, sigmaerr, wavelet, theta, nv, rpsigmaerr);
nsim = 500;
niter = 4;
vertcorr = 20;
horcorr = 25;
[Phimap, Claymap, Swmap, Time] = GeosPetroInversion3D(near, mid, far, TimeSeis, phiprior, clayprior, swprior, stdpetro, corrpetro, elastrain, petrotrain, wavelet, theta, nv, sigmaerr, vertcorr, horcorr, nsim, niter);
[X,Y,Z] = meshgrid((1:nil), (1:nxl), Time);


%% Plot results
figure(2)
subplot(131)
h = slice(X,Y,Z,Phimap,interx,intery,[]);
set(h,'edgecolor','none')
zlim([TimeSeis(1) TimeSeis(end)]); xlim([1 nil]); ylim([1 nxl]);
axis tight; grid on; box on; set(gca, 'ZDir', 'reverse')
xlabel('Inline'); ylabel('Crossline'); zlabel('Time (s)'); title('Porosity')
subplot(132)
h = slice(X,Y,Z,Claymap,interx,intery,[]);
set(h,'edgecolor','none')
zlim([TimeSeis(1) TimeSeis(end)]); xlim([1 nil]); ylim([1 nxl]);
axis tight; grid on; box on; set(gca, 'ZDir', 'reverse')
xlabel('Inline'); ylabel('Crossline'); zlabel('Time (s)'); title('Clay volume')
subplot(133)
h = slice(X,Y,Z,Swmap,interx,intery,[]);
set(h,'edgecolor','none')
zlim([TimeSeis(1) TimeSeis(end)]); xlim([1 nil]); ylim([1 nxl]);
axis tight; grid on; box on; set(gca, 'ZDir', 'reverse')
xlabel('Inline'); ylabel('Crossline'); zlabel('Time (s)'); title('Water saturation') 
