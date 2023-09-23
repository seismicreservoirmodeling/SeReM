close all; close all

%% Bayesian Petrophysical Inversion Driver %%
% In this script we apply the Bayesian petrophysical inversion method
% (Grana and Della Rossa, 2010) to predict the petrophysical properties 
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
varerr = 10^-4;
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
vppriormean = 4;
vspriormean = 2.4;
rhopriormean = 2.3;
Vpprior = vppriormean*ones(nxl,nil,nm);
Vsprior = vspriormean*ones(nxl,nil,nm);
Rhoprior = rhopriormean*ones(nxl,nil,nm);

%% Spatial correlation matrix
corrlength = 5*dt;
sigma0 = [ 0.0034    0.0037    0.0014
    0.0037    0.0042    0.0012
    0.0014    0.0012    0.0015];
sigmaprior = SpatialCovariance(corrlength, dt, nm, sigma0);

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
np = size(petrotrain,2);
elastrain = [VpTrain VsTrain RhoTrain];
nd = size(elastrain,2);
faciestrain = ones(size(PhiTrain));
faciestrain(PhiTrain>mean(PhiTrain(:)))=2;
nf = max(unique(faciestrain));

%% GMM
% grid petrophysical properties
ndiscr = 30;
phigrid = linspace(0.01, 0.4, ndiscr)';   
claygrid = linspace(0, 0.8, ndiscr)';   
swgrid = linspace(0, 1, ndiscr)'; 
% domain elastic properties
ndiscr = 25;
vpgrid  = linspace(3.2, 4.6, ndiscr)';
vsgrid = linspace(2, 3, ndiscr)';
rhogrid = linspace(2, 2.6, ndiscr)';
% Error
rpsigmaerr = 10^-2*eye(nd,nd);


%% Seismic inversion
[Vpmap, Vsmap, Rhomap, Phimap, Claymap, Swmap, Time] = BayesPetroInversion3D_GMM(near, mid, far, TimeSeis, Vpprior, Vsprior, Rhoprior, sigmaprior, elastrain, petrotrain, faciestrain, vpgrid, vsgrid, rhogrid, phigrid, claygrid, swgrid, sigmaerr, wavelet, theta, nv, rpsigmaerr);
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

