%% Rock physics inversion Driver %%
% In this script we apply the Bayesian Rock phyisics inversion to predict 
% the petrophysical properties (porosity, clay volume, and water saturation
% We implement 4 different options:
% Gaussian distribution and linear model 
% Gaussian mixture distribution and linear model (Grana, 2016)
% Gaussian mixture distribution and non-linear model (Grana and Della Rossa, 2010)
% Non-parametric distribution  and non-linear model (Grana, 2018).
% The linear rock physics model is a multi-linear regression and it is
% estimated from a training dataset.
% In this implementation of the non-linear model we assume that the joint
% distribution of model and % data can be estimated from a training dataset 
% (generated, for example, using a rock physics model)

%% Available data and parameters
% Load data (seismic data and time)
addpath(genpath('../SeReM'));
load Data/data4.mat

% training dataset
mtrain = [Phi Clay Sw];
nv = size(mtrain,2);
dtrain = [Vprpm Vsrpm Rhorpm];
nd = size(dtrain,2);
nf = max(unique(Facies));

% domain to evaluate the posterior PDF
phidomain = (0:0.005:0.4);   
cdomain = (0:0.01:0.8); 
swdomain = (0:0.01:1);
[P,V,S] = ndgrid(phidomain, cdomain, swdomain);
mdomain = [P(:) V(:) S(:)];

% measured data (elastic logs)
dcond = [Vp Vs Rho];
ns = size(dcond,1);

% matrix associated to the linear rock physics operator
R = zeros(nd,nv+1);
X = [mtrain ones(size(Phi))];
R(1,:) = regress(Vprpm,X); 
R(2,:) = regress(Vsrpm,X); 
R(3,:) = regress(Rhorpm,X); 

% Error
sigmaerr = 10^-2*eye(nd,nd);


%% Gaussian linear case
% prior model
mum = mean(mtrain);
sm = cov(mtrain);

% linearization
G = R(:,1:nv);
datacond = dcond-R(:,end)';

% inversion
[mupost, sigmapost, Ppost]  = RockPhysicsLinGaussInversion(mum, sm, G, mdomain, datacond, sigmaerr);

% posterior mean
Phipost = mupost(:,1);
Cpost = mupost(:,2);
Swpost = mupost(:,3);
Philp = mupost(:,1)-1.96*sqrt(sigmapost(1,1));
Clp = mupost(:,2)-1.96*sqrt(sigmapost(2,2));
Swlp = mupost(:,3)-1.96*sqrt(sigmapost(3,3));
Phiup = mupost(:,1)+1.96*sqrt(sigmapost(1,1));
Cup = mupost(:,2)+1.96*sqrt(sigmapost(2,2));
Swup = mupost(:,3)+1.96*sqrt(sigmapost(3,3));

% marginal posterior distributions
Ppostphi = zeros(ns,length(phidomain));
Ppostclay = zeros(ns,length(cdomain));
Ppostsw = zeros(ns,length(swdomain));
for i=1:ns
    Ppostjoint=reshape(Ppost(i,:),length(phidomain),length(cdomain),length(swdomain));
    Ppostphi(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),2);
    Ppostclay(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),3)),1);
    Ppostsw(i,:)=sum(squeeze(sum(squeeze(Ppostjoint),2)),1);
    Ppostphi(i,:)=Ppostphi(i,:)/sum(Ppostphi(i,:));
    Ppostclay(i,:)=Ppostclay(i,:)/sum(Ppostclay(i,:));
    Ppostsw(i,:)=Ppostsw(i,:)/sum(Ppostsw(i,:));
end

% plots
figure(1)
subplot(131)
pcolor(phidomain, Depth, Ppostphi); 
hold on; shading interp; colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
ylabel('Depth (m)'); xlabel('Porosity (v/v)');  xlim([0 0.4]);
plot(Phipost, Depth, 'r', 'LineWidth', 2);
subplot(132)
pcolor(cdomain, Depth, Ppostclay); 
hold on; shading interp; colorbar; 
plot(Clay, Depth, 'k', 'LineWidth', 2); 
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
plot(Cpost, Depth, 'r', 'LineWidth', 2);
subplot(133)
pcolor(swdomain, Depth, Ppostsw); 
hold on; shading interp; colorbar; 
plot(Sw, Depth, 'k', 'LineWidth', 2); 
plot(Swpost, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
hbc=colorbar; title(hbc, 'Probability');

%% Gaussian mixture linear case
% prior model
mum = zeros(nf,nv);
sm = zeros(nv,nv,nf);
pf = zeros(nf,1);
for k=1:nf
    pf(k) = sum(Facies==k)/ns;
    mum(k,:) = mean(mtrain(Facies==k,:));
    sm(:,:,k) = cov(mtrain(Facies==k,:));
end

[~, ~, ~, Ppost]  = RockPhysicsLinGaussMixInversion(pf, mum, sm, G, mdomain, datacond, sigmaerr);

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
figure(2)
subplot(131)
pcolor(phidomain, Depth, Ppostphi); 
hold on; shading interp; colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
ylabel('Depth (m)'); xlabel('Porosity (v/v)');  xlim([0 0.4]);
plot(Phimap, Depth, 'r', 'LineWidth', 2);
subplot(132)
pcolor(cdomain, Depth, Ppostclay); 
hold on; shading interp; colorbar; 
plot(Clay, Depth, 'k', 'LineWidth', 2); 
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
plot(Cmap, Depth, 'r', 'LineWidth', 2);
subplot(133)
pcolor(swdomain, Depth, Ppostsw); 
hold on; shading interp; colorbar; 
plot(Sw, Depth, 'k', 'LineWidth', 2); 
plot(Swmap, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
hbc=colorbar; title(hbc, 'Probability');


%% Gaussian mixture case
% The joint Gaussian mixture distribution is estimated from the training dataset
[~, ~, ~, Ppost] = RockPhysicsGaussMixInversion(Facies, mtrain, dtrain, mdomain, dcond, sigmaerr);
% The joint Gaussian distribution can also be used
% [mupost, sigmapost, Ppost] = RockPhysicsGaussInversion(mtrain, dtrain, mdomain, dcond, sigmaerr);

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
figure(3)
subplot(131)
pcolor(phidomain, Depth, Ppostphi); 
hold on; shading interp; colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
ylabel('Depth (m)'); xlabel('Porosity (v/v)');  xlim([0 0.4]);
plot(Phimap, Depth, 'r', 'LineWidth', 2);
subplot(132)
pcolor(cdomain, Depth, Ppostclay); 
hold on; shading interp; colorbar; 
plot(Clay, Depth, 'k', 'LineWidth', 2); 
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
plot(Cmap, Depth, 'r', 'LineWidth', 2);
subplot(133)
pcolor(swdomain, Depth, Ppostsw); 
hold on; shading interp; colorbar; 
plot(Sw, Depth, 'k', 'LineWidth', 2); 
plot(Swmap, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
hbc=colorbar; title(hbc, 'Probability');

%% Non-parametric case (Kernel density estimation)

% petrophysical domain discretization
ndiscr = 25;
phidomain = linspace(0, 0.4, ndiscr)';   
cdomain = linspace(0, 0.8, ndiscr)';   
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
figure(4)
subplot(131)
pcolor(phidomain, Depth, Ppostphi); 
hold on; shading interp; colorbar; 
plot(Phi, Depth, 'k', 'LineWidth', 2);  
ylabel('Depth (m)'); xlabel('Porosity (v/v)');  xlim([0 0.4]);
plot(Phimap, Depth, 'r', 'LineWidth', 2);
subplot(132)
pcolor(cdomain, Depth, Ppostclay); 
hold on; shading interp; colorbar; 
plot(Clay, Depth, 'k', 'LineWidth', 2); 
xlabel('Clay volume (v/v)'); xlim([0 0.8]);
plot(Cmap, Depth, 'r', 'LineWidth', 2);
subplot(133)
pcolor(swdomain, Depth, Ppostsw); 
hold on; shading interp; colorbar; 
plot(Sw, Depth, 'k', 'LineWidth', 2); 
plot(Swmap, Depth, 'r', 'LineWidth', 2);
xlabel('Water saturation (v/v)');  xlim([0 1]);
hbc=colorbar; title(hbc, 'Probability');
