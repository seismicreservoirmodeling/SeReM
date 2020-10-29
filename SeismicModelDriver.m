%% Seismic Model Driver %%
% In this script we apply the convolutional model of a wavelet and the
% linearized approximation of Zoeppritz equations to compute synthetic
% seismograms for different reflection angles
% The model parameterization is expressed in terms of P- and S-wave
% velocity and density.

%% Available data and parameters
% Load data (elastic properties and depth)
addpath(genpath('../SeReM/'))
load Data/data2.mat

%% Initial parameters
% number of variables
nv = 3;
% reflection angles 
theta = [15, 30, 45];

% travel time
dt = 0.001;
t0 = 1.8;
TimeLog = [t0; t0 + 2*cumsum(diff(Depth)./(Vp(2:end,:)))];
Time = (TimeLog(1):dt:TimeLog(end))'; 

% time-interpolated elastic log
Vp = interp1(TimeLog, Vp, Time);
Vs = interp1(TimeLog, Vs, Time);
Rho = interp1(TimeLog, Rho, Time);

% number of samples (seismic properties)
nd = length(Vp)-1;

%% Wavelet
% wavelet 
freq = 45;
ntw = 64;
[wavelet, tw] = RickerWavelet(freq, dt, ntw);


%% Plot elastic data
figure(1)
subplot(131)
plot(Vp, Time, 'k', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('P-wave velocity (km/s)'); ylabel('Time (s)');
subplot(132)
plot(Vs, Time, 'k', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('S-wave velocity (km/s)'); 
subplot(133)
plot(Rho, Time, 'k', 'LineWidth', 2);
axis tight; grid on; box on; set(gca, 'YDir', 'reverse');
xlabel('Density (g/cm^3)'); 

%% Synthetic seismic data
[Seis, TimeSeis] = SeismicModel (Vp, Vs, Rho, Time, theta, wavelet);
Snear = Seis(1:nd);
Smid = Seis(nd+1:2*nd);
Sfar = Seis(2*nd+1:end);


%% Plot seismic data
figure(2)
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