% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Takes polynomial weights and produces an M-Sequence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% coeffs (Px1 Integers) = Polynomial coefficients. For example, if the
% polynomial is D^5+D^3+D^1+1 then the coeffs vector will be [1;0;1;0;1;1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% MSeq (Wx1 Integers) = W bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mSeq]=fMSeqGen(coeffs)
degMax = length(coeffs) - 1;
nStates = 2 ^ (degMax) - 1;
regInit = ones(1, degMax);
regLast = regInit;
reg = zeros(nStates, degMax);
reg(1, :) = regInit;
temp=coeffs(2: end);
for iState = 2: nStates
reg(iState, 2: end) = regLast(1: end - 1);

reg(iState, 1) = mod(sum(regLast & temp), 2);
regLast = reg(iState, :);
end
mSeq = 1 - 2 * reg(:, end);
end