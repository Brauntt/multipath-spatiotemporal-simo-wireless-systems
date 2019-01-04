function [delayEst] = fChannelEstimation(symbolsOut, goldSeq)
% Function:
%   - perform channel estimation for the desired source using the received
%  signal
%
% InputArg(s):
%   - symbolsOut: channel symbol chips received
%   - goldSeq: gold sequence used in the modulation process
%
% OutputArg(s):
%   - delayEst: estimated delay of the signals
%
% Comments:
%   - this function is suitable for signals with only one path
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% obtain the maximum possible relative delay and number of signals
[nDelays, nSignals] = size(goldSeq);
% assume the max actual delay can be large (improve accuracy with cost of
% complexity)
delayMax = length(symbolsOut) - nDelays;
% estimated delay for all signals
delayEst = zeros(nSignals, 1);
% all possible correlation functions
corFun = zeros(delayMax, nSignals);
for iDelay = 1: delayMax
    % calculate the correlation functions for all possible delays
    corFun(iDelay, :) = abs(symbolsOut(iDelay: iDelay + nDelays - 1).' * goldSeq);
end
for iSignal = 1: nSignals
    % sort the correlation indexes
    [~, delayIndex] = sort(corFun(:, iSignal), 'descend');
    % find the residues of indexes and delete the repeated values
    uniqueResidue = unique(mod(delayIndex, nDelays), 'stable');
    % the first unique residue can suggest the delay
    delayEst(iSignal) = uniqueResidue(1) - 1;
end
% set the invalid estimation as zero
delayEst(delayEst < 0) = 0;
end
