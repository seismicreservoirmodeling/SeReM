function confmat = ConfusionMatrix(ftrue, fpred, nf)

% CONFUSION MATRIX computes the confusion matrix of a discrete
% classification
% INPUT ftrue = true model
%       fpred = predicted model
%       nf = number of possible outcomes (e.g. number of facies)
% OUTPUT confmat = confusion matrix (absolute frequencies)

% Written by Dario Grana (August, 2020)

ns = size(ftrue,1);
confmat = zeros(nf,nf);
for i=1:ns
    confmat(ftrue(i),fpred(i)) =  confmat(ftrue(i),fpred(i))+1;
end
