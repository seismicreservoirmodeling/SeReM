

%% Creat synthetic well data 
addpath(genpath('../SeReM/'))
WELL = read_las_file('Data\pseudowell_2856667402688.las');

Depth = WELL.curves(:,1);
facies = WELL.curves(:,2);
Phi = WELL.curves(:,3); Phi(Phi<=0.001) = 0.001;
v_clay = WELL.curves(:,5); v_clay(v_clay <=0.001) = 0.001;
sw = WELL.curves(:,4); sw (sw <=0.001) = 0.001;


[Vp, Vs, Rho] = RPM(Phi ,v_clay, sw); Vp = 1000*Vp; Vs = 1000*Vs;

std_vp = 100;
std_vs = 50;
std_rho = 0.05;
correlation_function = construct_correlation_function_beta(30,1,Vp,2);
Vp = Vp + std_vp*FFT_MA_3D(correlation_function, randn(size(Vp)));
Vs = Vs + std_vs*FFT_MA_3D(correlation_function, randn(size(Vp)));
Rho = Rho + std_rho*FFT_MA_3D(correlation_function, randn(size(Vp)));
AI = Vp.*Rho;

figure
subplot(121)
scatter(Phi, Vp, 25, facies);
subplot(122)
scatter(AI, Vp./Vs, 25, sw);

figure
subplot(121)
plot_histogram(AI, facies);
subplot(122)
plot_histogram(Vp./Vs, facies);



%% PRIOR SAMPLING
n_sim = 100;
mtrain = [];
dtrain = [];
for facie = 1:length(unique(facies))    
    idx_facie = (facies==facie);  
    Phi_train = mean(Phi(idx_facie)) + std(Phi(idx_facie))*randn(n_sim,1);
    v_clay_train = mean(v_clay(idx_facie)) + std(v_clay(idx_facie))*randn(n_sim,1);
    sw_train = mean(sw(idx_facie)) + std(sw(idx_facie))*randn(n_sim,1);        
    mtrain = [ mtrain; Phi_train, v_clay_train, sw_train];
    [Vp_train, Vs_train, Rho_train] = RPM(Phi_train, v_clay_train, sw_train); Vp_train = 1000*Vp_train; Vs_train = 1000*Vs_train;
    dtrain = [ dtrain; Vp_train, Vs_train, Rho_train];    
end



%% Non-parametric case (Kernel density estimation)

% petrophysical domain discretization
ndiscr = 25;
phidomain = linspace(0, 0.4, ndiscr)';   
cdomain = linspace(0, 1, ndiscr)';   
swdomain = linspace(0, 1, ndiscr)';   
mdomain = [phidomain cdomain swdomain];
% elastic domain discretization
vpdomain  = linspace(min(Vp), max(Vp),ndiscr)';
vsdomain = linspace(min(Vs), max(Vs),ndiscr)';
rhodomain = linspace(min(Rho), max(Rho),ndiscr)';
ddomain =[vpdomain vsdomain rhodomain];
% kernel bandwidths 
h = 5;
hm(1) = (max(phidomain)-min(phidomain))/h;
hm(2) = (max(cdomain)-min(cdomain))/h;
hm(3) = (max(swdomain)-min(swdomain))/h;
hd(1) = (max(vpdomain)-min(vpdomain))/h;
hd(2) = (max(vsdomain)-min(vsdomain))/h;
hd(3) = (max(rhodomain)-min(rhodomain))/h;


% measured data (elastic logs)
dcond = [Vp Vs Rho];
ns = size(dcond,1);

% inversion
Ppost = RockPhysicsKDEInversion(mtrain, dtrain, mdomain, ddomain, dcond, hm, hd);

% marginal posterior distributions
Ppostphi = zeros(ns,length(phidomain));
Ppostclay = zeros(ns,length(cdomain));
Ppostsw = zeros(ns,length(swdomain));
Phimap = zeros(ns,1);
Cmap = zeros(ns,1);
Swmap = zeros(ns,1);
for i=1:ns
    Ppostjoint=reshape(Ppost(i,:),length(phidomain),length(cdomain),length(swdomain));
    Ppostphi(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),2);
    Ppostclay(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),1);
    Ppostsw(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),2)),1);
    Ppostphi(i,:)=Ppostphi(i,:)/sum(Ppostphi(i,:));
    Ppostclay(i,:)=Ppostclay(i,:)/sum(Ppostclay(i,:));
    Ppostsw(i,:)=Ppostsw(i,:)/sum(Ppostsw(i,:));
    [~,Phimapind]=max(Ppostphi(i,:));
    [~,Cmapind]=max(Ppostclay(i,:));
    [~,Swmapind]=max(Ppostsw(i,:));
    Phimap(i)=phidomain(Phimapind);
    Cmap(i)=cdomain(Cmapind);
    Swmap(i)=swdomain(Swmapind);
end

% plots
figure
subplot(141)
imagesc(1, Depth, facies); 
ylabel('Depth (m)');
subplot(142)
pcolor(phidomain, Depth, Ppostphi); 
set(gca, 'YDir','reverse')
hold on; shading interp; %colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
xlabel('Porosity (v/v)');  xlim([0 0.4]);
plot(Phimap, Depth, 'r', 'LineWidth', 2);
subplot(143)
pcolor(cdomain, Depth, Ppostclay); 
set(gca, 'YDir','reverse')
hold on; shading interp; %colorbar; 
plot(v_clay, Depth, 'k', 'LineWidth', 2); 
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
plot(Cmap, Depth, 'r', 'LineWidth', 2);
subplot(144)
pcolor(swdomain, Depth, Ppostsw); 
set(gca, 'YDir','reverse')
hold on; shading interp; %colorbar; 
plot(sw, Depth, 'k', 'LineWidth', 2); 
plot(Swmap, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
