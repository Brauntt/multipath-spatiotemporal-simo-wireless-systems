% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform DS-QPSK Modulation on a vector of bits using a gold sequence
% with channel symbols set by a phase phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% bitsIn (Px1 Integers) = P bits of 1's and 0's to be modulated
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence to be used in the modulation process
% phi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% symbolsOut (Rx1 Complex) = R channel symbol chips after DS-QPSK Modulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [symbolsIn]=fDSQPSKModulator(bitsIn, goldSeq, phi)
% subStreams = 1 - 2 * (reshape(bitsIn, length(bitsIn) / 2, 2));
subStreams(:, 1) = 1 - 2 * bitsIn(1: 2: end);
subStreams(:, 2) = 1 - 2 * bitsIn(2: 2: end);
symbolsQpsk = 1 / sqrt(2) *(cos(phi) * subStreams(:, 1) + 1i * sin(phi) * subStreams(:, 2));
symbolsIn = symbolsQpsk * goldSeq';
symbolsIn = reshape(symbolsIn, numel(symbolsIn), 1);
end
