function [symbolsIn] = fDSQPSKModulator(bitsIn, goldSeq, phi)
% Function:
%   - perform DS-QPSK Modulation on a vector of bits using a gold sequence
%  with channel symbols set by a phase phi
%
% InputArg(s):
%   - bitsIn: bitsMax bits representing the image to be modulated
%   - goldSeq: gold sequence to be used in the modulation process
%   - phi: angle index in radian of the QPSK constellation points
%
% OutputArg(s):
%   - symbolsIn: channel symbol chips after DS-QPSK Modulation
%
% Comments:
%   - symbol power is set as 2 as required
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% ensure the bit is a integer times of 2
if mod(length(bitsIn), 2) ~= 0
    bitsIn(end + 1) = 0;
end
% number of symbols
nSymbols = length(bitsIn) / 2;
symbolsQpsk = zeros(nSymbols, 1);
% map bits pairs to symbols
for iSymbol = 1: nSymbols
    % retrieve bits pairs
    pair =  bitsIn(2 * iSymbol - 1: 2 * iSymbol);
    % map pairs to constellations
    if isequal(pair, [0; 0])
        symbolsQpsk(iSymbol) = sqrt(2) * (cos(phi) + 1i * sin(phi));
    elseif isequal(pair, [0; 1])
        symbolsQpsk(iSymbol) = sqrt(2) * (cos(phi + pi / 2) + 1i * sin(phi + pi / 2));
    elseif isequal(pair, [1; 1])
        symbolsQpsk(iSymbol) = sqrt(2) * (cos(phi + pi) + 1i * sin(phi + pi));
    else
        symbolsQpsk(iSymbol) = sqrt(2) * (cos(phi + 3 * pi / 2) + 1i * sin(phi + 3 * pi / 2));
    end
end
% encode the symbols by gold sequence
symbolsIn = symbolsQpsk * goldSeq';
% reshape the result symbols as a stream
symbolsIn = reshape(symbolsIn.', numel(symbolsIn), 1);
end
