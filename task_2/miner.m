function [shift] = miner(mSeq1, mSeq2, shiftMin)
shift = shiftMin - 1;
isBalanced = 0;
while(~isBalanced)
    shift = shift + 1;
    seq = 1 - 2 * (-mSeq1 == circshift(mSeq2, shift));
    isBalanced = sum(seq == -1) - sum(seq == 1) == 1; 
end
end

