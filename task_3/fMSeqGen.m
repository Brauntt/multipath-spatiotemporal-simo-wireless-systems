function [mSeq] = fMSeqGen(coeffs)
% Function:
%   - takes polynomial weights and produces an M-Sequence
%
% InputArg(s):
%   - coeffs: Polynomial coefficients. For example, if the polynomial is
%  D^5+D^3+D^1+1 then the coeffs vector will be [1;0;1;0;1;1]
%
% OutputArg(s):
%   - mSeq: M-Sequence of (-1)'s and 1's
%
% Comments:
%   - the input coefficients should be primitive for a valid result
%   - the result M-sequence is of (-1)'s and 1's
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% max degree of the polynomial
degMax = length(coeffs) - 1;
% number of possible states
nStates = 2 ^ (degMax) - 1;
% initial state of registers: the first is all one, others all zero
reg = zeros(nStates, degMax);
reg(1, :) = ones(1, degMax);
% state of the previous register
regPrev = reg(1, :);
% determine the following registers
for iState = 2: nStates
    % update the current register based on the previous and parity
    reg(iState, 2: end) = regPrev(1: end - 1);
    reg(iState, 1) = mod(sum(regPrev & coeffs(2: end)), 2);
    regPrev = reg(iState, :);
end
% map the result to (-1)'s and 1's
mSeq = 1 - 2 * reg(:, end);
end
