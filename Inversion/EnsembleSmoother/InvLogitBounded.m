function v = InvLogitBounded(w, minv, maxv)

% INVERSE LOGIT BOUNDED computes the inverse logit tranformation for
% bounded variables
% INPUT w = initial variable
%       minv = lower bound 
%       maxv = upper bound 
% OUTUPT index = transformed variable 

% Written by Dario Grana (August 2020)

% tranformation
v = (exp(w)*maxv+minv)./(exp(w)+1);
