function [doaEst, delayEst] = music(array, symbolsMatrix, covSample, goldSeq, nPaths)
% Function:
%   - find the direction of arrival based on MUSIC algorithm
%
% InputArg(s):
%   - array: coordinates of the receiving sensors
%   - covRx: covariance matrix of the received signal
%
% OutputArg(s):
%   - doa: direction of arrival in degree
%
% Comments:
%   - noise eigenvectors are orthogonal to array manifold
%   - the value of the cost function should approach zero at doa
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18
azimuth = 0: 180;
elevation = 0;
[nRelativeDelays, nSignals] = size(goldSeq);
nAnts = length(array);
nChips = length(goldSeq);
shiftMatrix = [zeros(1, 2 * nChips); eye(2 * nChips - 1) zeros(2 * nChips - 1, 1)];
goldSeqExtend = [goldSeq; zeros(size(goldSeq))];
costFunSub = zeros(length(azimuth), nRelativeDelays);
costFun = cell(nSignals, 1);
delayEst = cell(nSignals, 1);
doaEst = cell(nSignals, 1);
% Fourier transform matrix
phi = exp(-1i * pi / nChips);
phiVect = zeros(2 * nChips, 1);
ftMatrix = zeros(2 * nChips);
for iPos = 1: 2 * nChips
    phiVect(iPos) = phi ^ (iPos - 1);
end
for iPos = 1: 2 * nChips
    ftMatrix(:, iPos) = phiVect .^ (iPos - 1);
end
% z1
transformMatrix = kron(eye(nAnts), diag(ftMatrix * goldSeqExtend) \  ftMatrix);
% z1[n]
tfSignal = transformMatrix * symbolsMatrix;
nVectors = size(tfSignal,2);
covSmoothSignal = 0;
covSmoothTf = 0;
% Q = nPaths
subvectLength = 2 * nChips + 1 - nPaths;
phiSub = phiVect(1: subvectLength);
subVectSignal = cell(nAnts, nPaths);
subVectTf = cell(nAnts, nPaths);
vectSignal = reshape(tfSignal(:, 1), nAnts, 2 * nChips);
vectTf = reshape(diag(transformMatrix * transformMatrix'), nAnts, 2 * nChips);
for iAnt = 1: nAnts
    for iPath = 1: nPaths
        subVectSignal{iAnt, iPath} = vectSignal(iAnt, iPath: iPath + subvectLength - 1);
        subVectTf{iAnt, iPath} = vectTf(iAnt, iPath: iPath + subvectLength - 1);
    end
end
subVectSignal = cell2mat(subVectSignal');
subVectTf = cell2mat(subVectTf');
for iPath = 1: nPaths
    covSmoothSignal = covSmoothSignal + (1 / nPaths) * (subVectSignal(iPath, :)' * subVectSignal(iPath, :)) / length(subVectSignal(iPath, :));
    covSmoothTf = covSmoothTf + (1 / nPaths) * (subVectTf(iPath, :)' * subVectTf(iPath, :)) / length(subVectTf(iPath, :));
end
[nSources, eigVectSignal] = detection(covSmoothSignal);
% [nSources, eigVectTf] = detection(covSmoothTf);
covTf = transformMatrix * transformMatrix' / length(transformMatrix);
[nSources, eigVectTf] = detection(covTf);
for iSignal = 1: nSignals
    for iAzimuth = azimuth
        spvComponent = spv(array, [iAzimuth elevation]);
        for iDelay = 1: nRelativeDelays
            %         starManifold = kron(spvComponent, shiftMatrix ^ iDelay * goldSeqExtend);
            starManifold = kron(spvComponent, phiSub .^ iDelay);
            costFunSub(iAzimuth + 1, iDelay) = 1 ./ (starManifold' * fpoc([eigVectSignal eigVectTf]) * starManifold);
        end
    end
    costFun{iSignal} = costFunSub;
    [~, sortIndex] = sort(costFunSub(:), 'descend');
    tempIndex = sortIndex(1: nPaths(iSignal));
    [doaEst{iSignal}, delayEst{iSignal}] = ind2sub(size(costFunSub), tempIndex);
    doaEst{iSignal} = doaEst{iSignal} - 1;
end
end