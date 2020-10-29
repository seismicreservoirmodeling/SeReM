%% Rock Physics Model Driver %%
% In this script we apply different rock physics model to a synthetic well
% log of porosity. 
% The rock physics models include:
% Wyllie
% Raymer
% Soft sand model
% Stiff sand model
% Inclusion model for spherical pores
% Berryman inclusion model for ellptical pores
% We assume that the solid is 100% quartz and the fluid is 100% water. For
% mixtures of minerals and fluids, the elastic properties can be computed
% using the function Matrix Fluid Model.
% 

%% Available data and parameters
% Load data (porosity and depth)
addpath(genpath('../SeReM/'))
load Data/data1

% Initial parameters
Kmat=36; 
Gmat=45;
Rhomat=2.65;
Kfl=2.25;
Gfl=0;
Rhofl=1;


%% Empirical models
% Initial parameters
Vpmat = VelocityDefinitions(Kmat, Gmat, Rhomat);
Vpfl = VelocityDefinitions(Kfl, Gfl, Rhofl);

% Wyllie model
VpW = WyllieModel(Phi, Vpmat, Vpfl);
% Raymer model
VpR = RaymerModel(Phi, Vpmat, Vpfl);

% figures
figure(1)
plot(VpW, Depth, 'k', 'LineWidth', 2);
hold on
plot(VpR, Depth,'r', 'LineWidth', 2);
grid on; box on;
xlim([1.5 6.5]); xlabel('P-wave velocity (km/s)'); ylabel('Depth')
legend('Wyllie model', 'Raymer model')

figure(2)
scatter(Phi, VpW, 50, Phi, 'o');
hold on
scatter(Phi, VpR, 50, Phi, 'd');
grid on; box on;
xlim([0 0.3]); ylim([1.5 6.5]); 
xlabel('Porosity'); ylabel('P-wave velocity (km/s)');
legend('Wyllie model', 'Raymer model')

%% Granular media models
% Initial parameters
criticalporo=0.4;
coordnumber=7;
pressure=0.02;

% Density
Rho = DensityModel(Phi, Rhomat, Rhofl);

% Soft sand model
[VpSoft, VsSoft] = SoftsandModel(Phi, Rho, Kmat, Gmat, Kfl, criticalporo, coordnumber, pressure);
% Stiff sand model
[VpStiff, VsStiff] = StiffsandModel(Phi, Rho, Kmat, Gmat, Kfl, criticalporo, coordnumber, pressure);

% figures
figure(3)
subplot(121)
plot(VpSoft, Depth, 'k', 'LineWidth', 2);
hold on
plot(VpStiff, Depth,'r', 'LineWidth', 2);
grid on; box on;
xlim([1.5 6.5]); xlabel('P-wave velocity (km/s)'); ylabel('Depth')
subplot(122)
plot(VsSoft, Depth, 'k', 'LineWidth', 2);
hold on
plot(VsStiff, Depth, 'r', 'LineWidth', 2);
grid on; box on;
xlim([.5 4.5]); xlabel('S-wave velocity (km/s)'); ylabel('Depth');
legend('Soft sand model', 'Stiff sand model')

figure(4)
scatter(Phi, VpSoft, 50, Phi, 'o');
hold on
scatter(Phi, VpStiff, 50, Phi, 'd');
grid on; box on;
xlim([0 0.3]); ylim([1.5 6.5]); 
xlabel('Porosity'); ylabel('P-wave velocity (km/s)');
legend('Soft sand model', 'Stiff sand model')


%% Includion models 
% Initial parameters
Ar=0.2; % for elliptical inclusion model

% Density
Rho = DensityModel(Phi, Rhomat, Rhofl);

% Spherical inclusion model
[VpSph, VsSph] = SphericalInclusionModel(Phi, Rho, Kmat, Gmat, Kfl);
% Elliptical inclusion model
[VpEll, VsEll] = BerrymanInclusionModel(Phi, Rho, Kmat, Gmat, Kfl, Ar);

% figures
figure(5)
subplot(121)
plot(VpSph, Depth, 'k', 'LineWidth', 2);
hold on
plot(VpEll, Depth,'r', 'LineWidth', 2);
grid on; box on;
xlim([1.5 6.5]); xlabel('P-wave velocity (km/s)'); ylabel('Depth')
subplot(122)
plot(VsSph, Depth, 'k', 'LineWidth', 2);
hold on
plot(VsEll, Depth, 'r', 'LineWidth', 2);
grid on; box on;
xlim([.5 4.5]); xlabel('S-wave velocity (km/s)'); ylabel('Depth');
legend('Spherical pores', 'Elliptical pores (\alpha=0.2)')

figure(6)
scatter(Phi, VpSph, 50, Phi, 'o');
hold on
scatter(Phi, VpEll, 50, Phi, 'd');
grid on; box on;
xlim([0 0.3]); ylim([1.5 6.5]); 
xlabel('Porosity'); ylabel('P-wave velocity (km/s)');
legend('Spherical pores', 'Elliptical pores (\alpha=0.2)')