function l = RadialCorrLength(l, angles, theta, gamma)
%l = sqrt( (lmin^2*lmax^2) ./ (lmax^2*(sin(azim-theta)).^2 + lmin^2*(cos(azim-theta)).^2) );

% RADIAL CORR LENGTH computes the radial correlation length 
% INPUT l = 3x1 array of ranges for 3 dimensions
%       angles = 2x1 array os azimuth and dip angles
%       theta = matrix of radial coordinates
% OUTPUT l = radial correlation length

% Written by Dario Grana (August, 2020), updated by de Figueiredo (Sept, 2023)

azim = angles(1);
dip = angles(2);

if length(l)==2
    l = sqrt( (l(1)^2*l(2)^2) ./ (l(1)^2*(sin(azim-theta)).^2 + l(2)^2*(cos(azim-theta)).^2) );
else
    l = sqrt( (l(1)^2*l(2)^2*l(3)^2) ./ ( l(3)^2*(l(2)^2*(cos(azim-theta)).^2 + l(1)^2*(sin(azim-theta)).^2).*(cos(dip - gamma)).^2 + l(1)^2*l(2)^2*sin(dip - gamma).^2 ) );
end
   


