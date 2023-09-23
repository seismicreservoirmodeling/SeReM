function [Phimap, Claymap, Swmap, Time] = GeosPetroInversion3D(near, mid, far, TimeSeis, phiprior, clayprior, swprior, stdpetro, corrpetro, elastrain, petrotrain, wavelet, theta, nv, sigmaerr, vertcorr, horcorr, nsim, niter)

% GeosPetroInversion3D computes the posterior distribution of 
% petrophysical properties using the Ensemble Smoother method
% (Liu and Grana, 2018) with PFS simulations
% INPUT near = seismic near angle (nxl x nil x nd)
%       mid = seismic mid angle (nxl x nil x nd)
%       far = seismic far angle (nxl x nil x nd)
%       TimeSeis = seismic time vector (nd x 1)
%       phiprior = prior porosity model (nxl x nil x nm, with nm = nd+1)
%       clayprior = prior clay model (nxl x nil x nm, with nm = nd+1)
%       swprior = prior saturation model (nxl x nil x nm, with nm = nd+1)
%       stdpetro = vector with prior standard deviation of petrophysical 
%                   properties (3 x 1)
%       corrpetro = prior correlation matrix of petrophysical properties (3x3)
%       elastrain = matrix with training elastic data [Vp, Vs, density]
%                    (ns x 3)
%       petrotrain = matrix with training petrophysics data 
%                    [porosoty, clay, saturation](ns x 3)
%       wavelet = wavelet vector 
%       theta = vector of reflection angles 
%       nv = number of model variables
%       sigmaerr = covariance matrix of the error (nv*nsamples x nv*nsamples)
%       vertcorr = vertical correlation parameter
%       horcorr = horizontal correlation parameter
%       nsim = number of simulations
%       niter = number of ES assimiations 
% OUTPUT Phimap = Predicted porosity (nxl x nil x nm)
%       Claymap = Predicted clay (nxl x nil x nm)
%       Swmap = Predicted saturation (nxl x nil x nm)
%       Time =  time vector (nm x 1)

% Written by Dario Grana (June 2023)

nxl = size(near,1);
nil = size(near,2);
nd = size(near,3);
nm = nd+1;
dt = TimeSeis(2)-TimeSeis(1);
Time = (TimeSeis(1)-dt/2:dt:TimeSeis(end)+dt)';
Phimap = zeros(nxl, nil, nm);
Claymap = zeros(nxl, nil, nm);
Swmap = zeros(nxl, nil, nm);


%% Prior realizations
[Phisim, Claysim, Swsim] = ProbFieldSimulation3D(vertcorr,horcorr,phiprior,clayprior, swprior, stdpetro, corrpetro, nxl, nil, nm, nsim);

%% Rock physics model
regcoef = zeros(nv,nv+1);
regx = [petrotrain ones(size(petrotrain,1),1)];
regcoef(1,:) = regress(elastrain(:,1),regx); 
regcoef(2,:) = regress(elastrain(:,2),regx); 
regcoef(3,:) = regress(elastrain(:,3),regx); 


%% ESMDA petrophysical inversion
alpha = 1/niter;

for i=1:nxl
    disp(['Percentage progress of inversion: ', num2mstr(round(i/nxl*100)), ' %'])
    for j=1:nil
        PostModels = [squeeze(Phisim(i,j,:,:)); squeeze(Claysim(i,j,:,:)); squeeze(Swsim(i,j,:,:))];
        Phipost = PostModels(1:nm,:);
        Claypost = PostModels(nm+1:2*nm,:);
        Swpost = PostModels(2*nm+1:end,:);
        Phipost(Phipost<0)=0; Phipost(Phipost>0.4)=0.4;
        Claypost(Claypost<0)=0; Claypost(Claypost>0.8)=0.8;
        Swpost(Swpost<0)=0; Swpost(Swpost>1)=1;
        SeisData = [squeeze(near(i,j,:)); squeeze(mid(i,j,:)); squeeze(far(i,j,:))];
        SeisPred = ForwardGeopModel1D(Phipost, Claypost, Swpost, regcoef, Time, theta, wavelet, nm, nv, nsim);
        for h=1:niter
            [PostModels, ~] = EnsembleSmootherMDA(PostModels, SeisData, SeisPred, alpha, sigmaerr);
            Phipost = PostModels(1:nm,:);
            Claypost = PostModels(nm+1:2*nm,:);
            Swpost = PostModels(2*nm+1:end,:);
            Phipost(Phipost<0)=0; Phipost(Phipost>0.4)=0.4;
            Claypost(Claypost<0)=0; Claypost(Claypost>0.8)=0.8;
            Swpost(Swpost<0)=0; Swpost(Swpost>1)=1;
            SeisPred = ForwardGeopModel1D(Phipost, Claypost, Swpost, regcoef, Time, theta, wavelet, nm, nv, nsim);
        end
        % posterior models
        Phimap(i,j,:) = mean(Phipost,2);
        Claymap(i,j,:) = mean(Claypost,2);
        Swmap(i,j,:) = mean(Swpost,2);
    end
end


