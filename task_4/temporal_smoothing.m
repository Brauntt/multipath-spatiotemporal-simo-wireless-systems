function [covSmooth] = temporal_smoothing(nSubVects, lenSubVect, nAnts, nChips, obj)
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
%   - covSmooth: smoothed covariance matrix
%
% Comments:
%   - the diagram can be found in document
%
% Author & Date: Yang (i@snowztail.com) - 30 Dec 18

% number of samples
nSamples = size(obj, 2);
% subvector of certain sample
subVectPiece = cell(nSubVects, nAnts);
% all subvectors
subVect = cell(nSamples, nSubVects);
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
        subVect{iSample, iSubVect} = cell2mat(subVectPiece(iSubVect,:));
    end
end
for iSubVect = 1: nSubVects
    % convert elements to double
    temp = cell2mat(subVect(:, iSubVect));
    % calculate covariance matrices of subvectors
    covSubVect{iSubVect} = temp' * temp / length(temp);
end
% smoothed covariance matrix of object
covSmooth = mean(cat(3, covSubVect{:}), 3);
end

