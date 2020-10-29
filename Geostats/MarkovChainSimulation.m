function fsim = MarkovChainSimulation(T, ns, nsim)

% MARKOV CHAIN SIMULATION simulates 1D realizations of a discrete random 
% variable based on a stationary first-order Markov chain with given
% transition probability matrix 
% INPUT T = transition  probability matrix 
%       ns = number of samples
%       nsim = number of simulations
% OUTUPT fsim = realizations (ns, nsim)

% Written by Dario Grana (August 2020)

fsim = zeros(ns, nsim);
fprior=T^100; 
fprior=fprior(1,:);
for j=1:nsim
    fsim(1,j) = RandDisc(fprior);
    for i=2:ns
        fsim(i,j) = RandDisc(T(fsim(i-1,j),:));
    end
end
