function [nSources] = detector_mdl(nSamples, tfSignalSmooth)
% Function: 
%   - detector for practical sampled signal, based on minimum description
% length criterion
%
% InputArg(s):
%   - nSamples: number of samples
%   - tfSignalSmooth: smoothed signal
%
% OutputArg(s):
%   - nSources: number of sources
%
% Comments:
%   - source count equals minimum index of the function minus one
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18

[~, eigValue] = eig(tfSignalSmooth);
eigValue = sort(abs(diag(eigValue)));
nReceivers = length(eigValue);
mdlFun = (-nSamples) * (log(flip(cumprod(eigValue))) + (nReceivers: -1: 1)' .* (log((nReceivers: -1: 1)') - log(flip(cumsum(eigValue))))) + 1/2 * log(nReceivers) * (0: nReceivers - 1)' .* (2 * nReceivers: -1: nReceivers + 1)'; 
[~, minPos] = min(mdlFun);
nSources = minPos - 1;
end
