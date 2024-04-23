%% Seismic inversion Driver %%
% In this script we apply the Bayesian linearized AVO inversion method
% (Buland and Omre, 2003) to predict the elastic properties (P- and S-wave
% velocity and density) from seismic data.

%% Available data and parameters
% Load data (seismic data and time)
addpath(genpath('../SeReM/'))
load Data/data3.mat

%% Initial parameters
% number of samples (elastic properties)
nm = size(Snear,1)+1;
% reflection angles 
theta = [15, 30, 45];
ntheta = length(theta);
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

%% Spatial correlation matrix
% number of variables
%nv = 3; % For Aki & Richards
nv = 2; % For de Figueiredo & Grana
corrlength = 5*dt;
trow = repmat(0:dt:(nm-1)*dt,nm,1);
tcol = repmat((0:dt:(nm-1)*dt)',1,nm);
tdis = abs(trow-tcol);
sigmatime = exp(-(tdis./corrlength).^2);

%% Seismic inversion - Aki & Richards
nv = 3; % For Aki & Richards
sigma0 = cov([log(Vp),log(Vs),log(Rho)]);
sigmaprior = kron(sigma0, sigmatime);

Seis = [Snear; Smid; Sfar];
[mmap, mlp, mup, Time] = SeismicInversion(Seis, TimeSeis, Vpprior, Vsprior, Rhoprior, sigmaprior, sigmaerr, wavelet, theta, nv);

% Vpmap = mmap(1:nm);
% Vsmap = mmap(nm+1:2*nm);
% Rhomap = mmap(2*nm+1:end);
% Vplp = mlp(1:nm);
% Vslp = mlp(nm+1:2*nm);
% Rholp = mlp(2*nm+1:end);
% Vpup = mup(1:nm);
% Vsup = mup(nm+1:2*nm);
% Rhoup = mup(2*nm+1:end);
% 
% %% Plot results
% figure(2)
% subplot(131)
% plot(Vp, Time, 'k', 'LineWidth', 2);
% hold on;
% plot(Vpprior, Time, 'b', 'LineWidth', 1);
% plot(Vpmap, Time, 'r', 'LineWidth', 2);
% plot(Vplp, Time, 'r--', 'LineWidth', 1);
% plot(Vpup, Time, 'r--', 'LineWidth', 1);
% axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
% xlabel('P-wave velocity (km/s)'); ylabel('Time (s)');
% subplot(132)
% plot(Vs, Time, 'k', 'LineWidth', 2);
% hold on
% plot(Vsprior, Time, 'b', 'LineWidth', 1);
% plot(Vsmap, Time, 'r', 'LineWidth', 2);
% plot(Vslp, Time, 'r--', 'LineWidth', 1);
% plot(Vsup, Time, 'r--', 'LineWidth', 1);
% axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
% xlabel('S-wave velocity (km/s)'); 
% subplot(133)
% plot(Rho, Time, 'k', 'LineWidth', 2);
% hold on
% plot(Rhoprior, Time, 'b', 'LineWidth', 1);
% plot(Rhomap, Time, 'r', 'LineWidth', 2);
% plot(Rholp, Time, 'r--', 'LineWidth', 1);
% plot(Rhoup, Time, 'r--', 'LineWidth', 1);
% axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
% xlabel('Density (g/cm^3)'); 


%% Seismic inversion - de Figueiredo and Grana
nv = 2; % For de Figueiredo & Grana
sigma0 = cov([log(Rho.*Vp),log(Vp./Vs)]);
sigmaprior = kron(sigma0, sigmatime);

Seis = [Snear; Smid; Sfar];
[mmap, mlp, mup, Time] = SeismicInversion(Seis, TimeSeis, Vpprior, Vsprior, Rhoprior, sigmaprior, sigmaerr, wavelet, theta, nv);

Ipmap = mmap(1:nm);
VpVsmap = mmap(nm+1:2*nm);
Iplp = mlp(1:nm);
VpVslp = mlp(nm+1:2*nm);
Ipup = mup(1:nm);
VpVsup = mup(nm+1:2*nm);

%% Plot results
figure(3)
subplot(121)
plot(Vp.*Rho, Time, 'k', 'LineWidth', 2);
hold on;
plot(Vpprior.*Rhoprior, Time, 'b', 'LineWidth', 1);
plot(Ipmap, Time, 'r', 'LineWidth', 2);
plot(Iplp, Time, 'r--', 'LineWidth', 1);
plot(Ipup, Time, 'r--', 'LineWidth', 1);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('P-Impedance (km/s g/cm^3)'); ylabel('Time (s)');
subplot(122)
plot(Vp./Vs, Time, 'k', 'LineWidth', 2);
hold on
plot(Vpprior./Vsprior, Time, 'b', 'LineWidth', 1);
plot(VpVsmap, Time, 'r', 'LineWidth', 2);
plot(VpVslp, Time, 'r--', 'LineWidth', 1);
plot(VpVsup, Time, 'r--', 'LineWidth', 1);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Vp/Vs ratio'); 

