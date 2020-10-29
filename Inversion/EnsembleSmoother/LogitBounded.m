function w = LogitBounded(v, minv, maxv)

% LOGIT BOUNDED computes the logit tranformation for
% bounded variables
% INPUT v = initial variable
%       minv = lower bound of the domain
%       maxv = upper bound of the domain
% OUTUPT index = transformed variable 

% Written by Dario Grana (August 2020)
w = log((v-minv)./(maxv -v));
