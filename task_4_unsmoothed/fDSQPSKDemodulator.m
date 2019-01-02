function [bitsOut] = fDSQPSKDemodulator(symbolsOut, weightStRake, goldSeq, phi)
% Function:
%   - perform demodulation of the received data
%
% InputArg(s):
%   - symbolsOut: channel symbol chips received from each antenna
%   - weightStRake: the weight of the spatiotemporal rake beamformer
%   - goldSeq: gold sequence used in the modulation process
%   - phi: angle index in radian of the QPSK constellation points
%
% OutputArg(s):
%   - bitsOut: demodulated bits
%
% Comments:
%   - demap the symbols by angle rather than distance
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% combine the received streams on antennas by weights
symbolsOut = (weightStRake' * symbolsOut).';
% obtain the number of signals
nSignals = size(goldSeq, 2);
% symbol length
nSymbols = length(symbolsOut);
% output bits
bitsOut = zeros(2 * nSymbols, nSignals);
% symbol angles in QPSK modulation
angleQpsk = [phi, phi + pi / 2, phi - pi, phi - pi / 2];
% demodulate symbols by angle
for iSignal = 1: nSignals
    for iSymbol = 1: nSymbols
        % regard the symbol as the pattern with minimum angle difference
        [~, pos] = min(abs(angle(symbolsOut(iSymbol, iSignal)) - angleQpsk));
        if pos == 1
            bitsOut(2 * iSymbol - 1: 2 * iSymbol, iSignal) = [0; 0];
        elseif pos == 2
            bitsOut(2 * iSymbol - 1: 2 * iSymbol, iSignal) = [0; 1];
        elseif pos == 3
            bitsOut(2 * iSymbol - 1: 2 * iSymbol, iSignal) = [1; 1];
        else
            bitsOut(2 * iSymbol - 1: 2 * iSymbol, iSignal) = [1; 0];
        end
    end
end
end
