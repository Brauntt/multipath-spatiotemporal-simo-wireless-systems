function [doaEst, delayEst] = fChannelEstimation(array, symbolsOut, goldSeq, nPaths)
% Function:
%   - perform channel estimation for the desired source using the received
%  signal
%
% InputArg(s):
%   - array: array locations in half unit wavelength. If no array then
%  should be [0,0,0]
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
% obtain the maximum possible relative delay and number of signals
[nDelays, nSignals] = size(goldSeq);
% chip length
nChips = length(goldSeq);
% number of receiving antennas
nAnts = length(array);
% number of samples
nSamples = size(symbolsOut, 2);
% star manifold matrix
starMatrix = cell(nSignals);
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
% shifting matrix
shiftMatrix = [zeros(1, 2 * nChips); eye(2 * nChips - 1), zeros(2 * nChips - 1, 1)];
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
    % covariance matrix of transformed signal
    covTfSignal = tfSignal * tfSignal' / size(tfSignal, 2);
    % number of submatrices of spatial smoothing
    nSubMats = nPaths(iSignal);
    nSubMats = 2;
    % number of space-time subvectors
    nSubVects = nPaths(iSignal);
    % length of space-time subvectors (d + Q - 1 <= 2 * Nc)
    lenSubVect = 2 * nChips - 1 - nSubVects;
    % obtain Fourier transformation subvector
    ftSubVect = ftVect(1: lenSubVect);
    % then perform temporal smoothing of the spatial smoothed result
%     [covSmoothSignal] = temporal_smoothing(nSubVects, lenSubVect, nAnts, nChips, tfSignalSpatialSmooth);
%     [covSmoothSignal] = ts(tfSignalSpatialSmooth, nChips, nSubVects);

    % first perform spatial smoothing for signal
    [tfSignalSpatial] = spatial_smoothing(nSubMats, nAnts, covTfSignal);
    % then temporal smoothing
    [tfSignalSmooth] = temp(nSubVects, lenSubVect, nChips, tfSignalSpatial);
%     [tfSignalSmooth] = temporal_smoothing(nSubVects, lenSubVect, nAnts - nSubMats + 1, nChips, tfSignalSpatial);
    
    % first perform temporal smoothing for transformation
    tfMatrixTemporal = temp(nSubVects, lenSubVect, nChips, tfMatrix * tfMatrix');
%     [tfMatrixTemporal] = temporal_smoothing(nSubVects, lenSubVect, nAnts, nChips, tfMatrix * tfMatrix');
    % then spatial smoothing
    [tfMatrixSmooth] = spatial_smoothing(nSubMats, nAnts, tfMatrixTemporal);
%     [tfMatrixSpatial] = spatial_smoothing(nSubMats, nAnts, tfMatrix * tfMatrix');
%     [tfMatrixSmooth] = temporal_smoothing(nSubVects, lenSubVect, nAnts - nSubMats + 1, nChips, tfMatrixSpatial);
    
%     [covSmoothTf] = temporal_smoothing(nSubVects, lenSubVect, nAnts, nChips, tfMatrix * tfMatrix');
%     % smoothed covariance matrix of transformation
%     [covSmoothTf] = temporal_smoothing(nSubVects, lenSubVect, nAnts, nChips, tfMatrix * tfMatrix');
%     covTf = tfMatrix * tfMatrix' / length(tfMatrix);
%     nSub = length(covTf) - rank(covTf) + 1;
%     [tfSmooth] = spatial_smoothing(lenSubVect, nAnts, tfMatrix * tfMatrix');
    [nSourcesAic] = detector_aic(nSamples, tfSignalSmooth, tfMatrixSmooth);
    [doaEst{iSignal}, delayEst{iSignal}] = music(array(1: nAnts - nSubMats + 1, :), tfSignalSmooth, tfMatrixSmooth, nSourcesAic, ftSubVect, nDelays, nAnts, nSubMats, nPaths(iSignal));
end
% store the estimations in matrices as required
doaEst = [cell2mat(doaEst), zeros(length(cell2mat(doaEst)), 1)];
delayEst = cell2mat(delayEst);
% set the invalid estimation as zero
delayEst(delayEst < 0) = 0;
end
