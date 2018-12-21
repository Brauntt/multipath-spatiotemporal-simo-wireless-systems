function [bitsOut] = fDSQPSKDemodulator(symbolsOut, goldSeq, phi, delayEst, nPaths, fadingCoefs)
% Function:
%   - perform demodulation of the received data
%
% InputArg(s):
%   - symbolsOut: channel symbol chips received from each antenna
%   - goldSeq: gold sequence used in the modulation process
%   - phi: angle index in radian of the QPSK constellation points
%   - delayEst: estimated delay of the signals
%   - nPaths: number of paths for each source
%   - fadingCoefs: fading coefficients for each path
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
% weighted symbols by maximum ratio combining
symbolWeighted = zeros(nSymbols, nSignals);
% output bits
bitsOut = zeros(2 * nSymbols, nSignals);
% initialise the path counter
pathCounter = 1;
for iSignal = 1: nSignals
    % despreaded symbols
    symbolDespread = zeros(nSymbols, nPaths(iSignal));
    weightMrc = zeros(nPaths(iSignal), 1);
    for iPath = 1: nPaths(iSignal)
        % synced symbols
        symbolsSync = circshift(symbolsOut, -delayEst(pathCounter));
        % despread symbols
        symbolDespread(:, iPath) = reshape(symbolsSync(1: length(symbolsOut) - nDelays), nDelays, nSymbols).' * goldSeq(:, iSignal);
        % MRC weight is the complex conjugate of fading coefficients
        weightMrc(iPath) = fadingCoefs(pathCounter)';
        % update the path counter
        pathCounter = pathCounter + 1;
    end
    % calculate the weighted signal symbols
    symbolWeighted(:, iSignal) = sum(symbolDespread .* weightMrc.', 2);
end
% symbol angles in QPSK modulation
angleQpsk = [phi, phi + pi / 2, phi - pi, phi - pi / 2];
% demodulate symbols by angle
for iSignal = 1: nSignals
    for iSymbol = 1: nSymbols
        % regard the symbol as the pattern with minimum angle difference
        [~, pos] = min(abs(angle(symbolWeighted(iSymbol, iSignal)) - angleQpsk));
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
