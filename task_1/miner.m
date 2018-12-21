function [shift] = miner(mSeq1, mSeq2, shiftMin)
% Function:
%   - find the minimum shift that should be greater than a threshold for
%   two M-sequences to produce a balanced gold sequence
%
% InputArg(s):
%   - mSeq1, mSeq2: two M-sequences to produce gold sequences
%   - shiftMin: the desired minimum shift must be greater than or equal to
%   this threshold
%
% OutputArg(s):
%   - shift: the minimum shift of the second M-sequence to the right for a
%  balanced gold sequence
%
% Comments:
%   - this function only find the minimum shift and the balanced gold
%   sequence is not returned
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% ensure the check start from shiftMin
shift = shiftMin - 1;
% enable looping
isBalanced = 0;
% keep iterating until a balanced sequence is found
while (~isBalanced)
    shift = shift + 1;
    % produce gold sequence
    seq = 1 - 2 * (-mSeq1 == circshift(mSeq2, shift));
    % check whether balanced
    isBalanced = sum(seq == -1) - sum(seq == 1) == 1;
end
end

