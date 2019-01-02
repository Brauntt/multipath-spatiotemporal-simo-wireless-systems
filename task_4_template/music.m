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

% possible azimuth and elevation angles of arrival and delays
azimuth = 0: 180; elevation = 0; delay = 1: nDelays;
% cost function
costFun = zeros(length(azimuth), nDelays);
% obtain generalised noise eigenvectors
[eigVectNoise] = noise_detection(tfSignalSmooth, tfMatrixSmooth, nSources);
for iAzimuth = azimuth
    % the corresponding manifold vector
    spvComponent = spv(array, [iAzimuth elevation]);
    for iDelay = 1: nDelays
        % spatio-temporal array manifold
        starManifold = kron(spvComponent, ftSubVect .^ iDelay);
        % corresponding cost function
        costFun(iAzimuth + 1, iDelay) = 1 ./ (starManifold' * (eigVectNoise) * starManifold);
    end
end
% plot the 2-d MuSIC spectrum
plot2d3d(abs(costFun.'), azimuth, delay, 'dB', '2-Dimensional MuSIC Spectrum');
% find max value of cost function for all directions
[costFunDelay, doaIndex] = max(abs(costFun));
% then obtain positions of several max values that suggest delays
[~, delayEst] = maxk(costFunDelay, nPaths);
% sort by delay
delayEst = sort(delayEst(1: nPaths));
% and find corresponding doas
doaEst = doaIndex(delayEst) - 1;
end
