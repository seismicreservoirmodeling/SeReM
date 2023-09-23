function [Vpmap, Vsmap, Rhomap, Phimap, Claymap, Swmap, Time] = BayesPetroInversion3D_GMM(near, mid, far, TimeSeis, vpprior, vsprior, rhoprior, sigmaprior, elastrain, petrotrain, faciestrain, vpgrid, vsgrid, rhogrid, phigrid, claygrid, swgrid, sigmaerr, wavelet, theta, nv, rpsigmaerr)

% BayesPetroInversion3D_GMM computes the posterior distribution of 
% petrophysical properties using Bayesian linearized AVO inversion 
% (Buland and Omre, 2003), a Gaussian mixture rock physics likelihood
% (Grana and Della Rossa, 2010) and Chapman Kolmogorov equation
% INPUT near = seismic near angle (nxl x nil x nd)
%       mid = seismic mid angle (nxl x nil x nd)
%       far = seismic far angle (nxl x nil x nd)
%       TimeSeis = seismic time vector (nd x 1)
%       vpprior = prior Vp model (nxl x nil x nm, with nm = nd+1)
%       vsprior = prior Vs model (nxl x nil x nm, with nm = nd+1)
%       rhoprior = prior density model (nxl x nil x nm, with nm = nd+1)
%       elastrain = matrix with training elastic data [Vp, Vs, density]
%                    (ns x 3)
%       petrotrain = matrix with training petrophysics data 
%                    [porosoty, clay, saturation](ns x 3)
%       faciestrain = vector with training facies data 
%       vpgrid = vector of discretized Vp grid (ndiscr x 1)
%       vsgrid = vector of discretized Vs grid (ndiscr x 1)
%       rhogrid = vector of discretized density grid (ndiscr x 1)
%       phigrid = vector of discretized porosity grid (ndiscr x 1)
%       claygrid = vector of discretized clay grid (ndiscr x 1)
%       swgrid = vector of discretized satruration grid (ndiscr x 1)
%       sigmaerr = covariance matrix of the error (nv*nsamples x nv*nsamples)
%       wavelet = wavelet vector 
%       theta = vector of reflection angles 
%       nv = number of model variables
%       rpsigmaerr = rock physics  error variance
% OUTPUT Vpmap = Predicted Vp (nxl x nil x nm)
%       Vpmap = Predicted Vp (nxl x nil x nm)
%       Rhomap = Predicted density (nxl x nil x nm)
%       Phimap = Predicted porosity (nxl x nil x nm)
%       Claymap = Predicted clay (nxl x nil x nm)
%       Swmap = Predicted saturation (nxl x nil x nm)
%       Time =  time vector (nm x 1)

% Written by Dario Grana (June 2023)

nxl = size(near,1);
nil = size(near,2);
nm = size(near,3)+1;
ndiscr = length(phigrid);
dt = TimeSeis(2)-TimeSeis(1);
Time = (TimeSeis(1)-dt/2:dt:TimeSeis(end)+dt)';
Vpmap = zeros(nxl, nil, nm);
Vsmap = Vpmap;
Rhomap = Vpmap;
Phimap = Vpmap;
Claymap = Vpmap;
Swmap = Vpmap;

% rock physiscs likelihood
[pg,cg,sg] = ndgrid(phigrid, claygrid, swgrid);
petrogrid = [pg(:) cg(:) sg(:)];
[vpg,vsg,rg]=ndgrid(vpgrid, vsgrid, rhogrid);
elasgrid = [vpg(:) vsg(:) rg(:)];
Ppetro = RockPhysicsGMM(faciestrain, petrotrain, elastrain, petrogrid, elasgrid, rpsigmaerr);

% inversion
for i=1:nxl
    disp(['Percentage progress: ', num2mstr(round(i/nxl*100)), ' %'])
    for j=1:nil
        Seis = [reshape(near(i,j,:), nm-1, 1)  
            reshape(mid(i,j,:), nm-1, 1) 
            reshape(far(i,j,:), nm-1, 1)];
        [mmap, mtrans, strans, Time] = SeismicInversion3D(Seis, TimeSeis, squeeze(vpprior(i,j,:)), squeeze(vsprior(i,j,:)), squeeze(rhoprior(i,j,:)), sigmaprior, sigmaerr, wavelet, theta, nv);
        Vpmap(i,j,:) = mmap(1:nm);
        Vsmap(i,j,:) = mmap(nm+1:2*nm);
        Rhomap(i,j,:) = mmap(2*nm+1:end);
        vptrans = mtrans(1:nm);
        vstrans = mtrans(nm+1:2*nm);
        rhotrans = mtrans(2*nm+1:end);
        sigmatrans = [strans(round(nm/2),1) 0 0
                    0 strans(nm+round(nm/2),1) 0
                    0 0 strans(2*nm+round(nm/2),1)];
        Pseis = zeros(size(elasgrid,1), nm);
        for k=1:nm
            Pseis(:,k)=mvnpdf(elasgrid, [vptrans(k) vstrans(k) rhotrans(k)], sigmatrans);
            Pseis(:,k)=Pseis(:,k)/sum(Pseis(:,k));
        end
        Ppost = Ppetro' * Pseis;
        Ppostmarg = zeros(ndiscr,ndiscr,ndiscr,nm);
        Pphi = zeros(nm,ndiscr);
        Pclay = zeros(nm,ndiscr);
        Psw = zeros(nm,ndiscr);
        for k=1:nm
            Ppostmarg(:,:,:,k)=reshape(Ppost(:,k),ndiscr,ndiscr,ndiscr);
            Pphi(k,:)=sum(squeeze(sum(squeeze(Ppostmarg(:,:,:,k)),3)),2);
            Pclay(k,:)=sum(squeeze(sum(squeeze(Ppostmarg(:,:,:,k)),3)),1);
            Psw(k,:)=sum(squeeze(sum(squeeze(Ppostmarg(:,:,:,k)),2)),1);
            Pphi(k,:)=Pphi(k,:)/sum(Pphi(k,:));
            Pclay(k,:)=Pclay(k,:)/sum(Pclay(k,:));
            Psw(k,:)=Psw(k,:)/sum(Psw(k,:));
            [~,ii]=max(Pphi(k,:));
            [~,jj]=max(Pclay(k,:));
            [~,kk]=max(Psw(k,:));
            Phimap(i,j,k)=phigrid(ii);
            Claymap(i,j,k)=claygrid(jj);
            Swmap(i,j,k)=swgrid(kk);
        end
    end
end

