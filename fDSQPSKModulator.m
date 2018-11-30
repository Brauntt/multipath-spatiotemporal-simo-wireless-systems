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
% subStreams(:, 1) = 1 - 2 * bitsIn(1: 2: end);
% subStreams(:, 2) = 1 - 2 * bitsIn(2: 2: end);
% symbolsQpsk = 1 / sqrt(2) * (cos(phi) * subStreams(:, 1) + 1i * sin(phi) * subStreams(:, 2));
if mod(length(bitsIn), 2) ~= 0 
    bitsIn(end + 1) = 0;
end
nSymbols = length(bitsIn) / 2;
symbolsQpsk = zeros(nSymbols, 1);
for iSymbol = 1: nSymbols
   pair =  bitsIn(2 * iSymbol - 1: 2 * iSymbol);
   if isequal(pair, [0; 0])
       symbolsQpsk(iSymbol) = 1 / sqrt(2) * (cos(phi) + 1i * sin(phi));
   elseif isequal(pair, [0; 1])
       symbolsQpsk(iSymbol) = 1 / sqrt(2) * (cos(phi + pi / 2) + 1i * sin(phi + pi / 2));
   elseif isequal(pair, [1; 1])
       symbolsQpsk(iSymbol) = 1 / sqrt(2) * (cos(phi + pi) + 1i * sin(phi + pi));
   else
       symbolsQpsk(iSymbol) = 1 / sqrt(2) * (cos(phi + 3 * pi / 2) + 1i * sin(phi + 3 * pi / 2));
   end
end

symbolsIn = symbolsQpsk * goldSeq';
symbolsIn = reshape(symbolsIn.', numel(symbolsIn), 1);
end
