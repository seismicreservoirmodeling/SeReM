function sgsim = SeqGaussianSimulation(xcoords, dcoords, dvalues, xmean, xvar, l, type, krig, angles)

% SEQ GAUSSIAN SIMULATION  generates a realization of the random variable 
% conditioned on the available measurements using Sequential Gaussian
% Simulation
% INPUT xcoords = coordinates of the locations for the estimation (np, ndim)
%       dcoords = coordinates of measurements (nd,ndim)
%       dvalues = values of measurements (nd,1)
%       xmean = prior mean
%       xvar = prior variance
%       h = distance
%       correlation length, 1x1 for isotropic or 3x1 for anisotropic
%       type = function ype ('exp', 'gau', 'sph')
%       krig = kriging type (0=simple, 1=ordinary)
%       angles = angles for anisotropic variogram , 3x1  
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
    if krig == 0
        [krigmean, krigvar] = SimpleKriging(pathcoords(i,:), dc, dz, xmean, xvar, l, type, angles);
    else
        [krigmean, krigvar] = OrdinaryKriging(pathcoords(i,:), dc, dz, xvar, l, type, angles);
    end
    % realization
    if krigvar<0
        krigvar = 0;
    end
    simval(pathind(i)) = krigmean+sqrt(krigvar)*randn(1);
    % Adding simulated value the vector of conditioning data
    dcoords = [dcoords; pathcoords(i,:)];
    dvalues = [dvalues; simval(pathind(i))]; 
    %imagesc(simval)
    %drawnow
end
% Assigning the sampled values to the simulation grid
sgsim(sgsim==-999)=simval;
