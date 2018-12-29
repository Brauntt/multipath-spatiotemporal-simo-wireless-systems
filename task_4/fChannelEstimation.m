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

%% Initialisation
% possible azimuth and elevation angles of arrival
azimuth = 0: 180; elevation = 0;
% obtain the maximum possible relative delay and number of signals
[nDelays, nSignals] = size(goldSeq);
% chip length
nChips = length(goldSeq);
% number of receiving antennas
nAnts = length(array);
% star manifold matrix
starMatrix = cell(nSignals);
% cost function
costFun = zeros(length(azimuth), nDelays);
% Fourier transformation vector
ftVect = zeros(2 * nChips, 1);
% Fourier transformation matrix
ftMatrix = zeros(2 * nChips);
% estimated DOA of each signal
doaEst = cell(nSignals, 1);
% estimated delay of each signal
delayEst = cell(nSignals, 1);
% data vectorisation
[symbolsMatrix] = data_vectorisation(symbolsOut, nAnts, nChips);
% % number of samples
% nSamples = size(symbolsMatrix, 2);
% covariance matrix of symbol matrix
covSymbol = symbolsMatrix * symbolsMatrix' / length(symbolsMatrix);
% signal eigenvectors detection
[~, eigVectSignal] = detection(covSymbol);
% shifting matrix
shiftMatrix = [zeros(1, 2 * nChips); eye(2 * nChips - 1) zeros(2 * nChips - 1, 1)];
% extend the gold sequence by padding zeros to double length
goldSeqExtend = [goldSeq; zeros(size(goldSeq))];
% Fourier transformation constant
ftConst = exp(-1i * pi / nChips);
% construct Fourier transformation vector with ft constant
for iPos = 1: 2 * nChips
    ftVect(iPos) = ftConst ^ (iPos - 1);
end
% construct Fourier transformation matrix with ft vector
for iPos = 1: 2 * nChips
    ftMatrix(:, iPos) = ftVect .^ (iPos - 1);
end
%% delay and doa estimation
for iSignal = 1: nSignals
    % overall transformation matrix for certain signal
    tfMatrix = kron(eye(nAnts), diag(ftMatrix * goldSeqExtend(:, iSignal)) \  ftMatrix);
    % transformed signal
    tfSignal = tfMatrix * symbolsMatrix;
    % number of space-time subvectors
    nSubVects = nPaths(iSignal);    
    % length of space-time subvectors (d + Q - 1 <= 2 * Nc)
    lenSubVect = 2 * nChips + 1 - nSubVects;
    % smoothed covariance matrix of signal
    [covSmoothVect] = temporal_smooth(nSubVects, lenSubVect, nAnts, nChips, tfSignal);
    % smoothed covariance matrix of transformation
    [covSmoothTf] = temporal_smooth(nSubVects, lenSubVect, nAnts, nChips, tfMatrix * tfMatrix');
    for iAzimuth = azimuth
        % the corresponding manifold vector
        spvComponent = spv(array, [iAzimuth elevation]);
        for iDelay = 1: nDelays
            % spatio-temporal array manifold
            starManifold = kron(spvComponent, shiftMatrix ^ iDelay * goldSeqExtend(:, iSignal));
        end
    end
end
% store the estimations in matrices as required
doaEst = [cell2mat(doaEst), zeros(length(cell2mat(doaEst)), 1)];
delayEst = cell2mat(delayEst);
% set the invalid estimation as zero
delayEst(delayEst < 0) = 0;
end
