function [nSourcesEst, eigVectSignal] = detection(covMatrix, desiredNoisePower)
% Function: 
%   - detector for continous signal, based on eigendecomposition
%
% InputArg(s):
%   - covMatrix: covariance matrix of the received signal
%   - desiredNoisePower: noise power added to desired signal
%
% OutputArg(s):
%   - nSourcesEst: number of sources estimated
%   - eigVectorSignal: signal eigenvector
%
% Restraints:
%   - Traditional detection is based on observations, where the criterion
%   of regarding a eigenvalue as 'small' corresponding to noise is
%   determined by noise power. In this case we use a simplified model based
%   on comparison: the threshold is defined as double the noise power. The
%   precision of this function needs improvements.
%
% Comments:
%   - also obtain signal eigenvector to create subspace for MUSIC algorithm
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18

[eigVector, eigValue] = eig(covMatrix);
eigValue = abs(diag(eigValue));
% signal and noise eigenvalue threshold
eigNoiseThr = 2 * desiredNoisePower;
% estimated source number 
nSourcesEst = sum(eigValue > eigNoiseThr);
% signal eigenvector
eigVectSignal = eigVector(:, eigValue > eigNoiseThr);
end

