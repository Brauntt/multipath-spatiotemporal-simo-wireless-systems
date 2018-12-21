function [bitsOut] = fDSQPSKDemodulator(symbolsOut, goldSeq, phi, delayEst)
% Function:
%   - perform demodulation of the received data
%
% InputArg(s):
%   - symbolsOut: channel symbol chips received from each antenna
%   - goldSeq: gold sequence used in the modulation process
%   - phi: angle index in radian of the QPSK constellation points
%   - delayEst: estimated delay of the signals
%
% OutputArg(s):
%   - bitsOut: demodulated bits
%
% Comments:
%   - demap the symbols by angle rather than distance
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% obtain the maximum possible relative delay and number of signals
[nDelays, nSignals] = size(goldSeq);
% number of signal symbols
nSymbols = (length(symbolsOut) - nDelays) / nDelays;
% despreaded symbols
symbolDespread = zeros(nSymbols, nSignals);
% output bits
bitsOut = zeros(2 * nSymbols, nSignals);
for iSignal = 1: nSignals
    % synced symbols
    symbolsSync = circshift(symbolsOut, -delayEst(iSignal));
    % despread symbols with gold sequence
    symbolDespread(:, iSignal) = reshape(symbolsSync(1: length(symbolsOut) - nDelays), nDelays, nSymbols).' * goldSeq(:, iSignal);
end
% symbol angles in QPSK modulation
angleQpsk = [phi, phi + pi / 2, phi - pi, phi - pi / 2];
% demodulate symbols by angle
for iSignal = 1: nSignals
    for iSymbol = 1: nSymbols
        % regard the symbol as the pattern with minimum angle difference
        [~, pos] = min(abs(angle(symbolDespread(iSymbol, iSignal)) - angleQpsk));
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
