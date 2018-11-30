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

function [bitsOut]=fDSQPSKDemodulator(symbolsOut, goldSeq, phi)
% nDelays = length(goldSeq) = number of possible delays
[nDelays, nSignals] = size(goldSeq);
corFun = zeros(nDelays, nSignals);
symbolsDesp = zeros(length(symbolsOut) - nDelays, nSignals);
% temp = reshape(symbolsOut, numel(symbolsOut) / nDelays, nDelays);
for iDelay = 1: nDelays
    corFun(iDelay, :) = abs(symbolsOut(iDelay: iDelay + nDelays - 1).' * goldSeq);
%     corFun(iDelay) = abs(temp(iDelay, :) * goldSeq(:, 1));
end
[~, delay] = max(corFun);
delay = delay - 1;
nSymbols = (length(symbolsOut) - nDelays) / nDelays;
symbol = zeros(nSymbols, nSignals);
for iSignal = 1: nSignals
    temp = circshift(symbolsOut, -delay(iSignal));
    symbolsDesp(:, iSignal) = temp(1: length(symbolsOut) - nDelays);
%     symbol(:, iSignal) = reshape(symbolsDesp(:, iSignal), nSymbols, nDelays) * goldSeq(:, iSignal);
    symbol(:, iSignal) = reshape(symbolsDesp(:, iSignal), nDelays, nSymbols).' * goldSeq(:, iSignal);
end
bitsOut = zeros(2 * nSymbols, 1);
angleSymbols = [phi, phi + pi / 2, phi - pi, phi - pi / 2];
for iSymbol = 1: nSymbols
    [~, pos] = min(abs(angle(symbol(iSymbol, 1)) - angleSymbols));
    if pos == 1
    bitsOut(2 * iSymbol - 1: 2 * iSymbol) = [0; 0];
    elseif pos == 2
    bitsOut(2 * iSymbol - 1: 2 * iSymbol) = [0; 1];
    elseif pos == 3
    bitsOut(2 * iSymbol - 1: 2 * iSymbol) = [1; 1];
    else
    bitsOut(2 * iSymbol - 1: 2 * iSymbol) = [1; 0];
    end
end
% symbolsDesp1 = symbolsDesp(1: length(symbolsDesp) - max(delay(iSignal), iSignal));
% symbols1 = symbolsDesp1
% symbolsDesp = zeros(numel(symbolsOut)/length(goldSeq), nSignals);
% for iSignal = 1: nSignals
%     symbolsDesp(:, iSignal) = reshape(symbolsOut, numel(symbolsOut)/length(goldSeq), 15) * goldSeq(:, iSignal);
% end
flag = 1;
end
