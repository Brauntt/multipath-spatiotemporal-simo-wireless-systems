function [goldSeq] = fGoldSeq(mSeq1, mSeq2, shift)
% Function:
%   - takes two M-Sequences of the same length and produces a gold sequence
%  by adding a delay and performing modulo 2 addition
%
% InputArg(s):
%   - mSeq1, mSeq2: two M-sequences to produce gold sequences
%   - shift: number of chips to shift second M-Sequence to the right
%
% OutputArg(s):
%   - goldSeq: gold sequence of (-1)'s and 1's
%
% Comments:
%   - the result gold sequence is of (-1)'s and 1's
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% generate gold sequence by shifting the second M-sequence
goldSeq = 1 - 2 * (-mSeq1 == circshift(mSeq2, shift));
end
