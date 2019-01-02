function [objTemporalSmooth] = temporal_smoothing_diagram(nSubVects, lenSubVect, nAnts, nChips, obj)
% Function: 
%   - temporal smoothing
%
% InputArg(s):
%   - nSubVects: number of space-time subvectors
%   - lenSubVect: length of space-time subvectors (d + Q - 1 <= 2 * Nc)
%   - nAnts: number of receiving antennas
%   - nChips: chip length
%   - obj: object matrix to perform temporal smoothing
%
% OutputArg(s):
%   - objTemporalSmooth: temporally smoothed covariance matrix
%
% Comments:
%   - the diagram can be found in document
%
% Author & Date: Yang (i@snowztail.com) - 30 Dec 18

% number of samples
nSamples = size(obj, 2);
% subvector pieces
subVectPiece = cell(nSubVects, nAnts);
% subvectors set
subVectSet = cell(nSamples, nSubVects);
% all subvectors
subVect = cell(nSubVects, 1);
% covariance matrix of subvectors
covSubVect = cell(nSubVects, 1);
for iSample = 1: nSamples
    % transformed signal on each antenna
    tfSignalSplit = reshape(obj(:, iSample), nAnts, 2 * nChips);
    for iSubVect = 1: nSubVects
        for iAnt = 1: nAnts
            % obtain subvector piece
            subVectPiece{iSubVect, iAnt} = tfSignalSplit(iAnt, iSubVect: iSubVect + lenSubVect - 1);
        end
        % concatenate pieces for subvectors
        subVectSet{iSample, iSubVect} = cell2mat(subVectPiece(iSubVect,:));
        % convert elements to double
        subVect{iSubVect} = cell2mat(subVectSet(:, iSubVect));
    end
end
for iSubVect = 1: nSubVects
    % calculate covariance matrices of subvectors
    covSubVect{iSubVect} = subVect{iSubVect}' * subVect{iSubVect} / size(subVect{iSubVect}, 2);
end
% smoothed covariance matrix of object
objTemporalSmooth = mean(cat(3, covSubVect{:}), 3);
end

