% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performs channel estimation for the desired source using the received signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (Fx1 Complex) = R channel symbol chips received
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired source used in the modulation process!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% delay_estimate = Vector of estimates of the delays of each path of the
% desired signal
% DOA_estimate = Estimates of the azimuth and elevation of each path of the
% desired signal
% beta_estimate = Estimates of the fading coefficients of each path of the
% desired signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [delayEst] = fChannelEstimation(symbolsOut, goldSeq)
[nRelativeDelays, nSignals] = size(goldSeq);
delayMax = 10 * nRelativeDelays;
delayEst = zeros(nSignals, 1);
corFun = zeros(delayMax, nSignals);
for iDelay = 1: delayMax
    corFun(iDelay, :) = abs(symbolsOut(iDelay: iDelay + nRelativeDelays - 1).' * goldSeq);
end
for iSignal = 1: nSignals
% [~, delayIndex] = maxk(corFun(:, iSignal), nPaths(iSignal));
% delayEst(pathCounter: pathCounter + nPaths(iSignal) - 1) = sort(delayIndex) - 1;
[~, delayIndex] = sort(corFun(:, iSignal), 'descend');
temp = unique(mod(delayIndex, nRelativeDelays), 'stable');
delayEst(iSignal) = temp(1) - 1;
end
% [~, delayEst] = max(corFun);
% delayEst = delayEst - 1;
end
