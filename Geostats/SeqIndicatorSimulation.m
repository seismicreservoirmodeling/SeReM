function sgsim = SeqIndicatorSimulation(xcoords, dcoords, dvalues, nf, pprior, l, type)

% SEQ INDICATOR SIMULATION  generates a realization of the discrete random 
% variable conditioned on the available measurements using Sequential 
% Indicator Simulation
% INPUT xcoords = coordinates of the locations for the estimation (np, ndim)
%       dcoords = coordinates of measurements (nd,ndim)
%       dvalues = values of the measurements (ns, 1)
%       nf = number of possible outcomes (e.g. number of facies)
%       pprior = prior probability (1,nf)
%       h = distance
%       l = correlation length
%       type = function ype ('exp', 'gau', 'sph')
% OUTPUT sgsim = realization

% Written by Dario Grana (August, 2020)

% initial parameters
n = size(xcoords,1);
nd = size(dcoords,1);

% maximum number of conditioning data
nmax = 12;

% Data assignment to the simulation grid (-999 or measurements)
sgsim = -999*ones(n,1);
for i=1:nd
    [~,ind] = min(sum((xcoords-dcoords(i,:)).^2,2));
    sgsim(ind) = dvalues(i);
end

% random path of locations
np = n-nd;
nonsimcoords = xcoords(sgsim==-999,:);
pathind = randperm(np);
pathcoords = nonsimcoords(pathind,:);
simval = zeros(np,1);

% sequential simulation
for i=1:np
    if size(dcoords,1) < nmax
        dc = dcoords;
        dz  = dvalues;
    else
        % conditioning data selection 
        dv = sqrt(sum((dcoords-repmat(pathcoords(i,:),size(dcoords,1),1)).^2,2));
        [~,ind] = sort(dv);
        dc = dcoords(ind(1:nmax),:);
        dz = dvalues(ind(1:nmax));
    end
    % kriging
    [ikprob, ~] = IndicatorKriging(pathcoords(i,:), dc, dz, nf, pprior, l, type);
    % realization
    simval(pathind(i)) = RandDisc(ikprob);
    % Adding simulated value the vector of conditioning data
    dcoords = [dcoords; pathcoords(i,:)];
    dvalues = [dvalues; simval(pathind(i))]; 
end
% Assigning the sampled values to the simulation grid
sgsim(sgsim==-999)=simval;
