function Ppost  = RockPhysicsGMM(ftrain, mtrain, dtrain, mdomain, dcond, sigmaerr)

% RockPhysicsGMM computes the rock physics likelihood assuming
% Gaussian mixture distribution  from a training dataset 
% INPUT ftrain = vector with training facies data 
%       mtrain = matrix with training petrophysics data 
%                    [porosoty, clay, saturation](ns x 3)
%       dtrain = matrix with training elastic data [Vp, Vs, density]
%                    (ns x 3)
%       mdomain = petrophysical domain (created using ndgrid)
%       dcond = elastic domain (created using ndgrid)
%       rpsigmaerr = rock physics  error variance
% OUTUPT Ppost = petrophysical joint distribution 

% Written by Dario Grana (August 2020)
% Modified by Dario Grana (June 2023)

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
    sjoint(:,:,k) = cov(datatrain);
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
    if sum(pfpost(i,:))>0
        pfpost(i,:) = pfpost(i,:)/sum(pfpost(i,:));
    end
    lh=0;
    for k=1:nf
        lh=lh+pfpost(i,k)*mvnpdf(mdomain, mupost(i,:,k),sigmapost(:,:,k));
    end
    % posterior PDF
    if sum(lh>0)
        Ppost(i,:)=lh/sum(lh);
    else
        Ppost(i,:) = 0;
    end
end

