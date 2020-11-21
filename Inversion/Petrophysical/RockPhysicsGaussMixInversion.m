function [mupost, sigmapost, pfpost, Ppost]  = RockPhysicsGaussMixInversion(ftrain, mtrain, dtrain, mdomain, dcond, sigmaerr)

% ROCK PHYSICS GAUSS MIX INVERSION computes the posterior distribution of
% petrophysical properties conditioned on elastic properties assuming a
% Gaussian mixture distribution.
% The joint distribution of the Bayesian inversion approach is estimated
% from a training dataset 
% INPUT ftrain = training dataset of facies (ntrain, 1)
%       mtrain = training dataset of petrophysical properties (ntrain, nm)
%       dtrain = training dataset of elastic properties (ntrain, nd)
%       mdomain = discretized domain of petrophysical properties 
%                 (generated using meshgrid)
%       dcond = measured data (nsamples, nd)
%       sigmaerr = covariance matrix of the error (nd, nd)
% OUTUPT mupost = posterior mean (nsamples x nv, 1)
%        sigmapost = posterior covariance matrix (nv, nv) 
%        fpost = posterior weights (facies proportions) (nsamples, 1) 
%        Ppost = joint posterior distribution 

% initial parameters
nv = size(mtrain, 2);
nd = size(dtrain, 2);
nf = max(unique(ftrain));
ns = size(dcond,1);
datatrain = [mtrain dtrain];

% joint distribution
pf = zeros(nf,1);
mjoint = zeros(nf,nv+nd);
mum = zeros(nf,nv);
mud = zeros(nf,nd);
sjoint = zeros(nv+nd,nv+nd,nf);
sm = zeros(nv,nv,nf);
sd = zeros(nd,nd,nf);
smd = zeros(nv,nd,nf);
sdm = zeros(nd,nv,nf);
for k=1:nf
    pf(k) = sum(ftrain==k)/length(ftrain);
    mjoint(k,:) = mean(datatrain(ftrain==k,:));
    mum(k,:) = mjoint(k,1:nv);
    mud(k,:) = mjoint(k,nv+1:end);
    sjoint(:,:,k) = cov(datatrain(ftrain==k,:));
    sm(:,:,k) = sjoint(1:nv,1:nv,k);
    sd(:,:,k) = sjoint(nv+1:end,nv+1:end,k);
    smd(:,:,k) = sjoint(1:nv,nv+1:end,k);
    sdm(:,:,k) = sjoint(nv+1:end,1:nv,k);
end

% posterior distribution 
mupost = zeros(ns,nv,nf);
sigmapost = zeros(nv,nv,nf);
pfpost = zeros(ns,nf);
Ppost = zeros(ns,size(mdomain,1));
% posterior covariance matrices
for k=1:nf
    sigmapost(:,:,k) = sm(:,:,k)-smd(:,:,k)/(sd(:,:,k)+sigmaerr)*sdm(:,:,k);
%     [~,posdefcheck] = chol(sigmapost(:,:,k));
%     if posdefcheck~=0
%         [V,D]=eig(sigmapost(:,:,k));
%         d=diag(D);
%         d(d<=0)=eps;
%         sigmapost(:,:,k)= V*diag(d)*V';
%     end
end
for i=1:ns
    for k=1:nf
        % posterior means
        mupost(i,:,k) = mum(k,:)'+smd(:,:,k)/(sd(:,:,k)+sigmaerr)*(dcond(i,:)-mud(k,:))';
        % posterior weights
        pfpost(i,k) = pf(k)*mvnpdf(dcond(i,:),mud(k,:),sd(:,:,k));
    end
    den=sum(pfpost(i,:));
    lh=0;
    for k=1:nf
        pfpost(i,k)=pfpost(i,k)/den;
        lh=lh+pfpost(i,k)*mvnpdf(mdomain, mupost(i,:,k),sigmapost(:,:,k));
    end
    % posterior PDF
    Ppost(i,:)=lh/sum(lh);
end
