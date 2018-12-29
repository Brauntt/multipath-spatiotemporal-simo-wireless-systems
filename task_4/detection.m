function [nSources, eigVectSignal] = detection(covMatrix)
% Function: 
%   - detector for continous signal, based on eigendecomposition
%
% InputArg(s):
%   - covMatrix: covariance matrix of the received signal
%
% OutputArg(s):
%   - nSources: number of sources
%   - eigVectorSignal: signal eigenvector
%
% Restraints:
%   - Traditional detection is based on observations, where the criterion
%   of regarding a eigenvalue as 'small' corresponding to noise is
%   determined by human. In this case we use a simplified model that based
%   on comparison: the threshold is defined by the minimum eigenvalue and a
%   ratio. This ratio should be determined carefully. The precision of this
%   function needs improvements.
% Comments:
%   - also obtain signal eigenvector to create subspace for MUSIC algorithm
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18

[eigVector, eigValue] = eig(covMatrix);
eigValue = abs(diag(eigValue));
% assume max noise power / min noise power is below this ratio
% TO BE DESIGNED BY ACTUAL CASES
eigNoiseThr = 0.01;
nSources = sum(eigValue > eigNoiseThr);
eigVectSignal = eigVector(:, eigValue > eigNoiseThr); 
% do not use noise eigenvector directly for noise subspace
% noiseEigVector = eigVector(:, eigValue <= eigNoiseThr);
end

