function [objSpatialSmooth] = spatial_smoothing(nSubMats, nAnts, obj)
% Function:
%   - spatial smoothing
%
% InputArg(s):
%   - nSubMats: number of submatrices of spatial smoothing
%   - nAnts: number of receiving antennas
%   - obj: object matrix to perform spatial smoothing
%
% OutputArg(s):
%   - objSpatialSmooth: spatially smoothed covariance matrix
%
% Comments:
%   - the diagram can be found in document
%
% Author & Date: Yang (i@snowztail.com) - 30 Dec 18

% initialisation
objSpatialSmooth = 0;
% vector length on each antenna
lenVect = length(obj) / nAnts;
% size of smoothed object
lenSmooth = lenVect * (nAnts + 1 - nSubMats);
for iSubMat = 1: nSubMats
    % retrieve submatrices and sum them
    objSpatialSmooth = objSpatialSmooth + obj((iSubMat - 1) * lenVect + 1: (iSubMat - 1) * lenVect + lenSmooth, (iSubMat - 1) * lenVect + 1: (iSubMat - 1) * lenVect + lenSmooth);
end
% then take average
objSpatialSmooth = objSpatialSmooth / nSubMats;
end
