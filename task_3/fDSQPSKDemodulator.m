% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform demodulation of the received data using <INSERT TYPE OF RECEIVER>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% symbolsIn (Fx1 Integers) = R channel symbol chips received
% goldseq (Wx1 Integers) = W bits of 1's and 0's representing the gold
% sequence of the desired signal to be used in the demodulation process
% phi (Integer) = Angle index in degrees of the QPSK constellation points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% bitsOut (Px1 Integers) = P demodulated bits of 1's and 0's
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bitsOut] = fDSQPSKDemodulator(symbolsOut, weightSuperres, goldSeq, phi, delayEst, nPaths, fadingCoefs)
symbolsOut = symbolsOut * conj(weightSuperres);
% nRelativeDelays = length(goldSeq) = number of possible delays
[nRelativeDelays, nSignals] = size(goldSeq);
nSymbols = (length(symbolsOut) - nRelativeDelays) / nRelativeDelays;
symbol = zeros(nSymbols, nSignals);
pathCounter = 1;
for iSignal = 1: nSignals
    symbolCut = zeros(length(symbolsOut) - nRelativeDelays, nPaths(iSignal));
    symbolPath = zeros(nSymbols, nPaths(iSignal));
    weightMrc = zeros(nPaths(iSignal), 1);
    for iPath = 1: nPaths(iSignal)
        temp = circshift(symbolsOut, -delayEst(pathCounter));
        symbolCut(:, iPath) = temp(1: length(symbolsOut) - nRelativeDelays);
        symbolPath(:, iPath) = reshape(symbolCut(:, iPath), nRelativeDelays, nSymbols).' * goldSeq(:, iSignal);
        weightMrc(iPath) = fadingCoefs(pathCounter)';
%         weightMrc(iPath) = 1;
        pathCounter = pathCounter + 1;
    end
    symbol(:, iSignal) = sum(symbolPath .* weightMrc.', 2);
end
bitsOut = zeros(2 * nSymbols, nSignals);
angleSymbols = [phi, phi + pi / 2, phi - pi, phi - pi / 2];
for iSignal = 1: nSignals
    for iSymbol = 1: nSymbols
        [~, pos] = min(abs(angle(symbol(iSymbol, iSignal)) - angleSymbols));
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
