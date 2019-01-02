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


% possible azimuth and elevation angles of arrival and delays
azimuth = 0: 180; elevation = 0; delay = 1: nDelays;
% cost function
costFun = zeros(length(azimuth), nDelays);
% obtain generalised noise eigenvectors
[eigVectNoise] = noise_detection(tfSignalSmooth, tfMatrixSmooth, nSources);
% [eigVectNoise] = detect(tfSignalSmooth, nSources);
% spvComponent = spv(array(1: nAnts - nSubMats + 1, :), [azimuth' elevation * ones(length(azimuth), 1)]);
for iAzimuth = azimuth
    % the corresponding manifold vector
    spvComponent = spv(array, [iAzimuth elevation]);
% angles = [azimuth', ones(length(azimuth), 1) * elevation];
    for iDelay = 1: nDelays
        % spatio-temporal array manifold
        starManifold = kron(spvComponent, ftSubVect .^ iDelay);
%         -10*log10(real(diag(S'*EE*EE'*S)))'
%         costFun(:, iDelay) = -10 * log10(real(diag(starManifold' * (eigVectNoise * eigVectNoise') * starManifold)));
        costFun(iAzimuth + 1, iDelay) = 1 ./ starManifold' * (eigVectNoise * eigVectNoise') * starManifold;
%         costFun(:, iDelay) = 1 ./ (starManifold' * (eigVectNoise * eigVectNoise') * starManifold);
    end
end
% plot the 2-d MuSIC spectrum
plot2d3d(abs(costFun.'), azimuth, delay, 'dB', '2-Dimensional MuSIC Spectrum');
% find max value of cost function for all directions
[costFunDelay, doaIndex] = max(abs(costFun));
% then obtain positions of several max values that suggest delays
[~, delayEst] = maxk(costFunDelay, nPaths);
% sort by delay
delayEst = sort(delayEst);
% and find corresponding doas
doaEst = doaIndex(delayEst) - 1;
end
