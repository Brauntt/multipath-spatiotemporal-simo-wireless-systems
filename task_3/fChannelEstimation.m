function [doaEst, delayEst] = fChannelEstimation(array, symbolsOut, goldSeq, nPaths)
% Function:
%   - perform channel estimation for the desired source using the received
%  signal
%
% InputArg(s):
%   - symbolsOut: channel symbol chips received
%   - goldSeq: gold sequence used in the modulation process
%   - nPaths: number of paths for each source
%
% OutputArg(s):
%   - doaEst: estimated direction of arrival of the signal paths
%   - delayEst: estimated delay of the signal paths
%
% Comments:
%   - number of paths for each source should be known
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% possible azimuth and elevation angles of arrival
azimuth = 0: 180; elevation = 0;
% obtain the maximum possible relative delay and number of signals
[nDelays, nSignals] = size(goldSeq);
% chip length
nChips = length(goldSeq);
% number of receiving antennas
nAnts = length(array);
% cost function
costFun = zeros(length(azimuth), nDelays);
% estimated DOA of each signal
doaEst = cell(nSignals, 1);
% estimated delay of each signal
delayEst = cell(nSignals, 1);
% data vectorisation
[symbolsMatrix] = data_vectorisation(symbolsOut, nAnts, nChips);
% covariance matrix of symbol matrix
covSymbol = symbolsMatrix * symbolsMatrix' / length(symbolsMatrix);
% signal eigenvectors detection
[~, eigVectSignal] = detection(covSymbol);
% shifting matrix
shiftMatrix = [zeros(1, 2 * nChips); eye(2 * nChips - 1) zeros(2 * nChips - 1, 1)];
% extend the gold sequence by padding zeros to double length
goldSeqExtend = [goldSeq; zeros(size(goldSeq))];
for iSignal = 1: nSignals
    for iAzimuth = azimuth
        % the corresponding manifold vector
        spvComponent = spv(array, [iAzimuth elevation]);
        for iDelay = 1: nDelays
            % spatio-temporal array manifold
            starManifold = kron(spvComponent, shiftMatrix ^ iDelay * goldSeqExtend(:, iSignal));
            % the corresponding cost function
            costFun(iAzimuth + 1, iDelay) = 1 ./ (starManifold' * fpoc(eigVectSignal) * starManifold);
        end
    end
    % sort the cost function indexes
    [~, sortIndex] = sort(costFun(:), 'descend');
    % the few maximum 1-D indexes
    columnIndex = sortIndex(1: nPaths(iSignal));
    % convert to 2-D to obtain corresponding DOA and delays
    [doaEst{iSignal}, delayEst{iSignal}] = ind2sub(size(costFun), columnIndex);
    % convert indexes to real delays
    doaEst{iSignal} = doaEst{iSignal} - 1;
end
% store the estimations in matrices as required
doaEst = [cell2mat(doaEst), zeros(length(cell2mat(doaEst)), 1)];
delayEst = cell2mat(delayEst);
% set the invalid estimation as zero
delayEst(delayEst < 0) = 0;
end
