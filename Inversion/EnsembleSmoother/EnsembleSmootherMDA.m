function [PostModels, KalmanGain] = EnsembleSmootherMDA(PriorModels, SeisData, SeisPred, alpha, sigmaerr)

% ENSEMBLE SMOOTHER MDA computes the updated realizations of the
% model variables conditioned on the assimilated data using the 
% Ensemble Smoother Multiple Data Assimilation
% INPUT PriorModels = prior models realizations (nm, ne) 
%       SeisData = measured seismic data (nd, 1)
%       SeisPred = predicted data (nd, ne)
%       alpha = inflation coefficient 
%       sigmaerr = covariance matrix of the error (nd, nd)
% OUTPUT PostModels = updated models realizations (nm, ne) 
%        KalmanGain = Kalman Gain Matrix
    
% initial parameters
[nd, ne] = size(SeisPred);
% data perturbation
SeisPert = repmat(SeisData,1,ne) + sqrt(alpha*sigmaerr)*randn(nd, ne);
% mean models
mum = mean(PriorModels, 2);
mud = mean(SeisPred, 2);
% covariance matrices
smd = 1/(ne-1) * (PriorModels - mum) * (SeisPred - mud)';
sdd = 1/(ne-1) * (SeisPred - mud) * (SeisPred - mud)';
% Kalman Gain
KalmanGain = smd  * pinv(sdd + alpha*sigmaerr);
% Updated models
PostModels = PriorModels + KalmanGain*(SeisPert - SeisPred);