function [weightStRake] = spatiotemporal_rake(array, doaEst, delayEst, goldSeq, nPaths, fadingCoefs)
% Function: 
%   - obtain the directional gains of the spatiotemporal rake beamformer
%
% InputArg(s):
%   - array: coordinates of the receiving sensors
%   - doaEst: directions of the desired user
%   - delayEst: delays of the desired user
%   - goldSeq: gold sequence used in the modulation process
%   - nPaths: number of paths for each source
%   - fadingCoefs: path fading coefficients
%
% OutputArg(s):
%   - weightStRake: the weight on receiving antennas
%
% Comments:
%   - spatiotemporal rake beamformer is suitable for single user system
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18

% number of signals (should be 1 for single user system) and number of
% chips
[nChips, nSignals] = size(goldSeq);
% number of antennas
nAnts = length(array);
% extend the gold sequence by padding zeros to double length
goldSeqExtend = [goldSeq; zeros(size(goldSeq))];
% shifting matrix
shiftMatrix = [zeros(1, 2 * nChips); eye(2 * nChips - 1) zeros(2 * nChips - 1, 1)];
% path counter
pathCounter = 0;
for iSignal = 1: nSignals
    % spatio-temporal array manifold
    starManifold = zeros(2 * nAnts * nChips, nPaths(iSignal));
    for iPath = 1: nPaths(iSignal)
        % update path counter
        pathCounter = pathCounter + 1;
        % array manifold of certain path (doa)
        spvComponent = spv(array, doaEst(pathCounter, :));
        % obtain STAR manifold by gold sequence and array manifold
        starManifold(:, iPath) = kron(spvComponent, shiftMatrix ^ delayEst(pathCounter) * goldSeqExtend(:, iSignal));
    end
end
% weights of the spatiotemporal rake beamformer
weightStRake = starManifold * fadingCoefs;
end

