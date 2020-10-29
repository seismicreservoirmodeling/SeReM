%% Facies Classification Driver %%
% In this script we illustrate Bayesian facies classification using 
% two assumptions: 
% Example 1: Gaussian model
% Example 2: non-parametric model estimated using KDE

% available data
addpath(genpath('../SeReM/'))
load Data/data4.mat
data = [Vp Rho];
ns = size(data,1);
nv = size(data,2);

% domain
v = (3:0.005:5);
r = (2:0.01:2.8);
[V, R] = meshgrid(v,r);
domain = [V(:) R(:)];

%% Gaussian model (2 components)
nf = max(unique(Facies));
fp = zeros(nf,1);
mup = zeros(nf,nv);
sp = zeros(nv,nv,nf);
for k=1:nf
    fp(k) = sum(Facies==k)/ns;
    mup(k,:) = mean(data(Facies==k,:));
    sp(:,:,k) = cov(data(Facies==k,:));
end

% likelihood function
GaussLikeFun = zeros(length(r),length(v),nf);
for k=1:nf
    lf = mvnpdf(domain,mup(k,:),sp(:,:,k)); 
    lf = lf/sum(lf);
    GaussLikeFun(:,:,k) = reshape(lf,length(r),length(v));
end

% plot likelihood
figure(1)
plot(data(:,1), data(:,2), '.k')
hold on
for k=1:nf
    contour(V,R,GaussLikeFun(:,:,k), 'LineWidth', 2);
end
grid on; box on;
xlabel('P-wave velocity (km/s)'); ylabel('Density (g/cm^3)');

% classification
[fmap, fpost] = BayesGaussFaciesClass(data, fp, mup, sp);

% confusion matrix (absolute frequencies)
confmat = ConfusionMatrix(Facies, fmap, nf);
% reconstruction rate 
reconstrate = confmat./repmat(sum(confmat,2),1,nf);
% recognition rate 
recognrate = confmat./repmat(sum(confmat,1),nf,1);
% estimation index  
estimindex = reconstrate-recognrate;  %#ok<*NASGU>

% plot results
figure(2)
subplot(141)
plot(Vp, Depth, 'k', 'LineWidth', 2);
xlabel('P-wave velocity (km/s)'); ylabel('Depth (m)');
grid on; box on; set(gca, 'YDir', 'reverse')
subplot(142)
plot(Rho, Depth, 'k', 'LineWidth', 2);
xlabel('Density (g/cm^3)'); set(gca, 'YTickLabel', []);
grid on; box on; set(gca, 'YDir', 'reverse')
subplot(143)
plot(fpost,Depth, 'LineWidth', 2);
xlabel('Facies probability'); set(gca, 'YTickLabel', []); legend('Sand', 'Shale');
grid on; box on; set(gca, 'YDir', 'reverse')
subplot(144)
imagesc([],Depth,fmap);
xlabel('Predicted facies'); set(gca, 'YTickLabel', []);
grid on; box on;

%% Non parametric model (2 components)
% training data
dtrain = data;
ftrain = Facies;

% likelihood function
KDElikefun = zeros(length(r),length(v),nf);
for k=1:nf
    lf = ksdensity(data(Facies==k,:),domain, 'Kernel', 'epanechnikov'); 
    lf = lf/sum(lf);
    KDElikefun(:,:,k) = reshape(lf,length(r),length(v));
end

% plot likelihood
figure(3)
plot(data(:,1), data(:,2), '.k')
hold on
for k=1:nf
    contour(V,R,KDElikefun(:,:,k), 'LineWidth', 2);
end
grid on; box on;
xlabel('P-wave velocity (km/s)'); ylabel('Density (g/cm^3)');

% classification
[fmap, fpost] = BayesKDEFaciesClass(data, dtrain, ftrain, fp, domain);

% confusion matrix (absolute frequencies)
confmat = ConfusionMatrix(Facies, fmap, nf);
% reconstruction rate 
reconstrate = confmat./repmat(sum(confmat,2),1,nf);
% recognition rate 
recognrate = confmat./repmat(sum(confmat,1),nf,1);
% estimation index  
estimindex = reconstrate-recognrate;

% plot results
figure(4)
subplot(141)
plot(Vp, Depth, 'k', 'LineWidth', 2);
xlabel('P-wave velocity (km/s)'); ylabel('Depth (m)');
grid on; box on; set(gca, 'YDir', 'reverse')
subplot(142)
plot(Rho, Depth, 'k', 'LineWidth', 2);
xlabel('Density (g/cm^3)'); set(gca, 'YTickLabel', []);
grid on; box on; set(gca, 'YDir', 'reverse')
subplot(143)
plot(fpost,Depth, 'LineWidth', 2);
xlabel('Facies probability'); set(gca, 'YTickLabel', []); legend('Sand', 'Shale');
grid on; box on; set(gca, 'YDir', 'reverse')
subplot(144)
imagesc([],Depth,fmap);
xlabel('Predicted facies'); set(gca, 'YTickLabel', []);
grid on; box on;

