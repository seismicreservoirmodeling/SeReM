function [mupost, sigmapost, pfpost, Ppost]  = RockPhysicsLinGaussMixInversion(pf, mum, sm, G, mdomain, dcond, sigmaerr)

% ROCK PHYSICS LINEAR GAUSS MIX INVERSION computes the posterior 
% distribution petrophysical properties conditioned on elastic properties 
% assuming a Gaussian mixture distribution and a linear rock physics model
% INPUT pf =  prior weights (facies proportions) (nf, 1)
%       mum = prior means of petrophysical propertiees (nf, nv)
%       sm = prior covariance matrices of petrophysical propertiees (nv, nv, nf)
%       G = rock physics operator matrix
%       mdomain = discretized domain of petrophysical properties 
%                 (generated using meshgrid)
%       dcond = measured data (nsamples, nd)
%       sigmaerr = covariance matrix of the error (nd, nd)
% OUTUPT mupost = posterior mean (nsamples, nv, nf)
%        sigmapost = posterior covariance matrix (nv, nv, nf) 
%        fpost = posterior weights (nsamples, nf) 
%        Ppost = joint posterior distribution 

% Written by Dario Grana (August 2020)

% initial parameters
nv = size(mum,2);
nf = size(mum,1);
nd = size(dcond, 2);
ns = size(dcond,1);

% analytical calculations
mud = zeros(nf,nd);
sd = zeros(nd,nd,nf);
smd = zeros(nv,nd,nf);
sdm = zeros(nd,nv,nf);
for k=1:nf
    mud(k,:) = G*mum(k,:)';
    sd(:,:,k) = G*sm(:,:,k)*G';
    smd(:,:,k) = sm(:,:,k)*G';
    sdm(:,:,k) = G*sm(:,:,k);
end

% posterior distribution
mupost = zeros(ns,nv,nf);
sigmapost = zeros(nv,nv,nf);
pfpost = zeros(ns,nf);
Ppost = zeros(ns,size(mdomain,1));
% posterior covariance matrices
for k=1:nf
    sigmapost(:,:,k) = sm(:,:,k)-smd(:,:,k)/(sd(:,:,k)+sigmaerr)*sdm(:,:,k);
end
% analytical solution
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
