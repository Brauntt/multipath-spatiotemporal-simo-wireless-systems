function [doaEst, delayEst] = music(array, covRx, goldSeq, nPaths)
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
nChips = length(goldSeq);
costFun = cell(nSignals, 1);
costFunSub = zeros(length(azimuth), nRelativeDelays);
delayEst = cell(nSignals, 1);
doaEst = cell(nSignals, 1);
% costFun = zeros(nSignals, length(azimuth), nRelativeDelays);
% costFun = zeros(length(azimuth), 1);
% starManifold = cell(nRelativeDelays, 1);
[~, eigVectSignal] = detection(covRx);
shiftMatrix = [zeros(1, 2 * nChips); eye(2 * nChips - 1) zeros(2 * nChips - 1, 1)];
goldSeqExtend = [goldSeq; zeros(size(goldSeq))];
% goldSeqExtend = repelem(goldSeqExtend, 1, nPaths');
for iSignal = 1: nSignals
    for iAzimuth = azimuth
        spvComponent = spv(array, [iAzimuth elevation]);
        %     starManifold = kron(spvComponent, shiftMatrix) * goldSeqExtend;
        %     for iSignal = 1: nSignals
        %     costFun(iSignal, iAzimuthAngle + 1) = 1 ./ (starManifold(:, iSignal)' * fpoc(eigVectSignal) * starManifold(:, iSignal));
        %     end
        for iDelay = 1: nRelativeDelays
%             starManifold = kron(spvComponent, shiftMatrix ^ iDelay) * goldSeqExtend(:, iSignal);
            starManifold = kron(spvComponent, shiftMatrix ^ iDelay * goldSeqExtend(:, iSignal));
            costFunSub(iAzimuth + 1, iDelay) = 1 ./ (starManifold' * fpoc(eigVectSignal) * starManifold);
        end
    end
    costFun{iSignal} = costFunSub;
    [~, sortIndex] = sort(costFunSub(:), 'descend');
    tempIndex = sortIndex(1: nPaths(iSignal));
    [doaEst{iSignal}, delayEst{iSignal}] = ind2sub(size(costFunSub), tempIndex);
    doaEst{iSignal} = doaEst{iSignal} - 1;
end
doaEst = [cell2mat(doaEst), zeros(length(cell2mat(doaEst)), 1)];
delayEst = cell2mat(delayEst);
end
