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

function [delayEst] = fChannelEstimation(symbolsOut, goldSeq, nPaths)
[nDelays, nSignals] = size(goldSeq);
corFun = zeros(nDelays, nSignals);
for iDelay = 1: nDelays
    corFun(iDelay, :) = abs(symbolsOut(iDelay: iDelay + nDelays - 1).' * goldSeq);
end
pathCounter = 1;
delayEst = zeros(sum(nPaths), 1);
for iSignal = 1: nSignals
[~, delayIndex] = maxk(corFun(:, iSignal), nPaths(iSignal));
delayEst(pathCounter: pathCounter + nPaths(iSignal) - 1) = sort(delayIndex) - 1;
pathCounter = pathCounter + nPaths(iSignal);
% delaySet{iSignal} = sort(delayIndex) - 1;
end
end

% for iSignal = 1: nSignals
%    for iPath = 1: nPaths(iSignal)
%        corFun = zeros(nDelays, 1);
%        for iDelay = 1: nDelays
%            corFun(iDelay) = abs(symbolsOut(iDelay: iDelay + nDelays - 1).' * goldSeq(:, iSignal));
%        end
%        [~, delayEst] = maxk(corFun, );
%    end
% end
% end
