function index = RandDisc(p)

% RANDDISC samples a discrete random variable with a given probability
% mass function
% INPUT p = probabilities 
% OUTUPT index = sampled value

% Written by Dario Grana (August 2020)

u=rand(1);
index=1;
s=p(1);
while ((u>s) && (index<length(p)))
    index=index+1;
    s=s+p(index);
end
