%% Ensemble Smoother Seismic Inversion Driver %%
% In this script we apply the Ensemble Smoother inversion method
% (Liu and Grana, 2018) to predict the elastic properties (P- and S-wave
% velocity and density) from seismic data.

%% Available data and parameters
% Load data (seismic data, reference elastic properties, and time)
addpath(genpath('../SeReM/'))
load Data/data3.mat

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
Vpprior = filtfilt(b, a, Vp);
Vsprior = filtfilt(b, a, Vs);
Rhoprior = filtfilt(b, a, Rho);
mprior = [Vpprior Vsprior Rhoprior];

%% Spatial correlation matrix
corrlength = 5*dt;
trow = repmat(0:dt:(nm-1)*dt,nm,1);
tcol = repmat((0:dt:(nm-1)*dt)',1,nm);
tdis = trow-tcol;
sigmatime = exp(-(tdis./corrlength).^2);
sigma0 = cov([Vp,Vs,Rho]);


%% Prior realizations
Vpsim = zeros(nm, nsim);
Vssim = zeros(nm, nsim);
Rhosim = zeros(nm, nsim);
SeisPred = zeros(nd*ntheta, nsim);
for i=1:nsim  
    msim = CorrelatedSimulation(mprior, sigma0, sigmatime);
    Vpsim(:,i) = msim(:, 1);
    Vssim(:,i) = msim(:, 2);
    Rhosim(:,i) = msim(:, 3);  
    [SeisPred(:,i), TimeSeis] = SeismicModel (Vpsim(:,i), Vssim(:,i), Rhosim(:,i), Time, theta, wavelet);
end

% plot of prior models
figure(2)
subplot(131)
plot(Vpsim, Time, 'b', 'LineWidth', 1);
hold on;
plot(Vp, Time, 'k', 'LineWidth', 2);
plot(Vpprior, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('P-wave velocity (km/s)'); ylabel('Time (s)');
subplot(132)
plot(Vssim, Time, 'b', 'LineWidth', 1);
hold on;
plot(Vs, Time, 'k', 'LineWidth', 2);
plot(Vsprior, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('S-wave velocity (km/s)'); 
subplot(133)
plot(Rhosim, Time, 'b', 'LineWidth', 1);
hold on;
plot(Rho, Time, 'k', 'LineWidth', 2);
plot(Rhoprior, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Density (g/cm^3)'); 


%% ESMDA seismic inversion
niter = 4;
alpha = 1/niter;  % sum alpha = 1
PriorModels = [Vpsim; Vssim; Rhosim];
SeisData = [Snear; Smid; Sfar];
PostModels = PriorModels;
for j=1:niter
    [PostModels, KalmanGain] = EnsembleSmootherMDA(PostModels, SeisData, SeisPred, alpha, sigmaerr);
    Vppost = PostModels(1:nm,:);
    Vspost = PostModels(nm+1:2*nm,:);
    Rhopost = PostModels(2*nm+1:end,:); 
    for i=1:nsim 
        [SeisPred(:,i), TimeSeis] = SeismicModel (Vppost(:,i), Vspost(:,i), Rhopost(:,i), Time, theta, wavelet);
    end
end

% posterior mean models
mpost = mean(PostModels, 2);
Vpmean = mpost(1:nm);
Vsmean = mpost(nm+1:2*nm);
Rhomean = mpost(2*nm+1:end);

%% Plot results
figure(3)
subplot(131)
plot(Vppost, Time, 'b', 'LineWidth', 1);
hold on;
plot(Vp, Time, 'k', 'LineWidth', 2);
plot(Vpmean, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('P-wave velocity (km/s)'); ylabel('Time (s)');
subplot(132)
plot(Vspost, Time, 'b', 'LineWidth', 1);
hold on;
plot(Vs, Time, 'k', 'LineWidth', 2);
plot(Vsmean, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('S-wave velocity (km/s)'); 
subplot(133)
plot(Rhopost, Time, 'b', 'LineWidth', 1);
hold on;
plot(Rho, Time, 'k', 'LineWidth', 2);
plot(Rhomean, Time, 'r', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Density (g/cm^3)'); 
