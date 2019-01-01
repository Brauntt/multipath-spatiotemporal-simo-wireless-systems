function [doaEst, delayEst] = music(array, tfSignalSmooth, tfMatrixSmooth, nSources, ftSubVect, nDelays, nAnts, nSubMats, nPaths)
% Function:
%   - find the direction of arrival based on MUSIC algorithm
%
% InputArg(s):
%   - array: coordinates of the receiving sensors
%   - covMatrix: covariance matrix of the received signal
%
% OutputArg(s):
%   - doa: direction of arrival in degree
%
% Comments:
%   - noise eigenvectors are orthogonal to array manifold
%   - the value of the cost function should approach zero at doa
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18

% azimuth = 0: 180;
% elevation = 0;
% costFun = zeros(length(azimuth), 1);
% [nSourses, eigVectorSignal] = detection(covMatrix);
% for iAzimuthAngle = azimuth
%     spvComponent = spv(array, [iAzimuthAngle elevation], mainlobe);
%     costFun(iAzimuthAngle + 1) = spvComponent' * fpoc(eigVectorSignal) * spvComponent;
% end
% [~, doaIndex] = mink(costFun, nSourses);
% doaAzimuth = doaIndex - 1;
% doa = [doaAzimuth elevation * ones(size(doaAzimuth))];
% doa = sortrows(doa, 1);


% possible azimuth and elevation angles of arrival
azimuth = 0: 180; elevation = 0;
% cost function
costFun = zeros(length(azimuth), nDelays);
% obtain generalised noise eigenvectors
[eigVectNoise] = detection(tfSignalSmooth, tfMatrixSmooth, nSources);
% spvComponent = spv(array(1: nAnts - nSubMats + 1, :), [azimuth' elevation * ones(length(azimuth), 1)]);
for iAzimuth = azimuth
    % the corresponding manifold vector
    spvComponent = spv(array(1: nAnts - nSubMats + 1, :), [iAzimuth elevation]);
    for iDelay = 1: nDelays
        % spatio-temporal array manifold
        %             starManifold = kron(spvComponent, shiftMatrix ^ iDelay * goldSeqExtend(:, iSignal));
        starManifold = kron(spvComponent, ftSubVect .^ iDelay);
        costFun(iAzimuth + 1, iDelay) = 1 ./ (starManifold' * (eigVectNoise * eigVectNoise') * starManifold);
%         costFun(:, iDelay) = 1 ./ (starManifold' * (eigVectNoise * eigVectNoise') * starManifold);
    end
end
% sort the cost function indexes
[~, sortIndex] = sort(costFun(:), 'descend');
% the few maximum 1-D indexes
columnIndex = sortIndex(1: nPaths);
% convert to 2-D to obtain corresponding DOA and delays
[doaEst, delayEst] = ind2sub(size(costFun), columnIndex);
% convert indexes to real delays
doaEst = doaEst - 1;
end
