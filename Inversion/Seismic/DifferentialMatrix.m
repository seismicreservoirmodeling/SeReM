function D=DifferentialMatrix(nt,nv)

% DIFFERENTIAL MATRIX computes the differential matrix for discrete
% differentiation
% INPUT nt = number of samples
%       nv = numbr of model variables
% OUTUPT D = differential matrix

% Written by Dario Grana (August 2020)

I  = eye(nt);
B  = zeros(nt,nt);
B(2:end,1:end-1) = -eye(nt-1);
I  = (I+B); I=I(2:end,:);
for i=1:nv
    D((i-1)*(nt-1)+1:i*(nt-1),(i-1)*nt+1:i*nt)=I;
end
