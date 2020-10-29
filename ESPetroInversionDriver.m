%% Ensemble Smoother Petrophysical Inversion Driver %%
% In this script we apply the Ensemble Smoother inversion method
% (Liu and Grana, 2018) to predict the petrophysical properties 
% from seismic data.

%% Available data and parameters
% Load data (seismic data, reference petroelastic properties, and time)
addpath(genpath('../SeReM/'))
load Data/data5.mat

%% Initial parameters
% number of samples (elastic properties)
nm = size(Snear,1)+1;
% number of samples (seismic data)
nd = size(Snear,1);
% number of variables
nv = 3;
% reflection angles 
theta = [15, 30, 45];
ntheta = length(theta);
% time sampling
dt = TimeSeis(2)-TimeSeis(1);
% error variance
varerr = 10^-4;
sigmaerr = varerr*eye(ntheta*(nm-1));
% number of realizations
nsim = 500;

%% Wavelet
% wavelet 
freq = 45;
ntw = 64;
[wavelet, tw] = RickerWavelet(freq, dt, ntw);

%% Linear Rock physics model
R = zeros(nd,nv+1);
X = [Phi Clay Sw ones(size(Phi))];
R(1,:) = regress(Vp,X); 
R(2,:) = regress(Vs,X); 
R(3,:) = regress(Rho,X); 

%% Plot seismic data
figure(1)
subplot(131)
plot(Snear, TimeSeis, 'k', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Near'); ylabel('Time (s)');
subplot(132)
plot(Smid, TimeSeis, 'k', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Mid'); 
subplot(133)
plot(Sfar, TimeSeis, 'k', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Far'); 

%% Prior model (filtered well logs)
nfilt = 3;
cutofffr = 0.04;
[b, a] = butter(nfilt, cutofffr);
Phiprior = filtfilt(b, a, Phi);
Clayprior = filtfilt(b, a, Clay);
Swprior = filtfilt(b, a, Sw);
mprior = [Phiprior Clayprior Swprior];

%% Spatial correlation matrix
corrlength = 5*dt;
trow = repmat(0:dt:(nm-1)*dt,nm,1);
tcol = repmat((0:dt:(nm-1)*dt)',1,nm);
tdis = trow-tcol;
sigmatime = exp(-(tdis./corrlength).^2);
sigma0 = cov([Phi,Clay,Sw]);


%% Prior realizations
Phisim = zeros(nm, nsim);
Claysim = zeros(nm, nsim);
Swsim = zeros(nm, nsim);
SeisPred = zeros(nd*ntheta, nsim);
for i=1:nsim   
    msim = CorrelatedSimulation(mprior, sigma0, sigmatime);
    Phisim(:,i) = msim(:,1);
    Claysim(:,i) = msim(:,2);
    Swsim(:,i) = msim(:,3);
end
Phisim(Phisim<0)=0; Phisim(Phisim>0.4)=0.4;
Claysim(Claysim<0)=0; Claysim(Claysim>0.8)=0.8;
Swsim(Swsim<0)=0; Swsim(Swsim>1)=1;
[Vpsim, Vssim, Rhosim] = LinearizedRockPhysicsModel(Phisim, Claysim, Swsim, R);
for i=1:nsim   
    [SeisPred(:,i), TimeSeis] = SeismicModel (Vpsim(:,i), Vssim(:,i), Rhosim(:,i), Time, theta, wavelet);
end

% plot of prior models
figure(2)
subplot(131)
plot(Phisim, Time, 'b', 'LineWidth', 1);
hold on;
plot(Phi, Time, 'k', 'LineWidth', 2);
plot(Phiprior, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Porosity'); ylabel('Time (s)');
subplot(132)
plot(Claysim, Time, 'b', 'LineWidth', 1);
hold on;
plot(Clay, Time, 'k', 'LineWidth', 2);
plot(Clayprior, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Clay volume'); 
subplot(133)
plot(Swsim, Time, 'b', 'LineWidth', 1);
hold on;
plot(Sw, Time, 'k', 'LineWidth', 2);
plot(Swprior, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Water saturation'); 


%% ESMDA petrophysical inversion
niter = 4;
alpha = 1/niter;
PriorModels = [Phisim; Claysim; Swsim];
SeisData = [Snear; Smid; Sfar];
SeisPred = zeros(nd*ntheta, nsim);
PostModels = PriorModels;
for j=1:niter
    [PostModels, KalmanGain] = EnsembleSmootherMDA(PostModels, SeisData, SeisPred, alpha, sigmaerr);
    Phipost = PostModels(1:nm,:); 
    Claypost = PostModels(nm+1:2*nm,:);
    Swpost = PostModels(2*nm+1:end,:); 
    Phipost(Phipost<0)=0; Phipost(Phipost>0.4)=0.4;
    Claypost(Claypost<0)=0; Claypost(Claypost>0.8)=0.8;
    Swpost(Swpost<0)=0; Swpost(Swpost>1)=1;
    [Vppost, Vspost, Rhopost] = LinearizedRockPhysicsModel(Phipost, Claypost, Swpost, R);
    for i=1:nsim 
        [SeisPred(:,i), TimeSeis] = SeismicModel (Vppost(:,i), Vspost(:,i), Rhopost(:,i), Time, theta, wavelet);
    end
end

% posterior mean models
mpost = mean(PostModels, 2);
Phimean = mpost(1:nm);
Claymean = mpost(nm+1:2*nm);
Swmean = mpost(2*nm+1:end);

%% Plot results
figure(3)
subplot(131)
plot(Phipost, Time, 'b', 'LineWidth', 1);
hold on;
plot(Phi, Time, 'k', 'LineWidth', 2);
plot(Phimean, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Porosity'); ylabel('Time (s)');
subplot(132)
plot(Claypost, Time, 'b', 'LineWidth', 1);
hold on;
plot(Clay, Time, 'k', 'LineWidth', 2);
plot(Claymean, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Clay volume'); 
subplot(133)
plot(Swpost, Time, 'b', 'LineWidth', 1);
hold on;
plot(Sw, Time, 'k', 'LineWidth', 2);
plot(Swmean, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Water saturation'); 
