function [objTemporalSmooth] = temporal_smoothing(nSubVects, lenSubVect, nChips, obj)
% Function:
%   - temporal smoothing of covariance matrix of sampled signal
%
% InputArg(s):
%   - nSubVects: number of space-time subvectors
%   - lenSubVect: length of space-time subvectors (d + Q - 1 <= 2 * Nc)
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

% number of segments
nSegs = length(obj) / (2 * nChips);
for iSegRow = 1: nSegs
    for iSegCol = 1: nSegs
        % data on each antenna
        objSeg = obj((iSegRow - 1) * 2 * nChips + 1: iSegRow * 2 * nChips, (iSegCol - 1) * 2 * nChips + 1: iSegCol * 2 * nChips);
        objSub = 0;
        for iSubVect = 1: nSubVects
            % divide into subvectors and sum them
            objSub = objSub + objSeg(iSubVect: iSubVect + lenSubVect - 1, iSubVect: iSubVect + lenSubVect - 1);
        end
        % take average of subvectors
        objSub = objSub / nSubVects;
        % the submatrix is part of smoothed result
        objTemporalSmooth((iSegRow - 1) * lenSubVect + 1: iSegRow * lenSubVect, (iSegCol - 1) * lenSubVect + 1: iSegCol * lenSubVect) = objSub;
    end
end
end
