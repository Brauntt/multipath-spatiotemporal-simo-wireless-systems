function [covSmooth] = temporal_smooth(nSubVects, lenSubVect, nAnts, nChips, object)

% number of samples
nSamples = size(object, 2);
% subvector of certain sample
subVectPiece = cell(nSubVects, nAnts);
% all subvectors
subVect = cell(nSamples, nSubVects);
% covariance matrix of subvectors
covSubVect = cell(nSubVects, 1);
for iSample = 1: nSamples
    % transformed signal on each antenna
    tfSignalSplit = reshape(object(:, iSample), nAnts, 2 * nChips);
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
    % calculate covariance matrices of subvectors
    covSubVect{iSubVect} = cov(cell2mat(subVect(:, iSubVect)));
end
% smoothed covariance matrix of object
covSmooth = mean(cat(3, covSubVect{:}), 3);
end

