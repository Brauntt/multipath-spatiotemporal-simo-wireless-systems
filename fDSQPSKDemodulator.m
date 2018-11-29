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
nDelays = 100;
[temp, nSignals] = size(goldSeq);
for iDelay = 1: nDelays
    fun(iDelay) = abs(symbolsOut(iDelay:iDelay+temp-1).'*goldSeq(:,1));
end
[~, delay] = max(fun);
de = mod(delay, 15)
% symbolsDesp = zeros(numel(symbolsOut)/length(goldSeq), nSignals);
% for iSignal = 1: nSignals
%     symbolsDesp(:, iSignal) = reshape(symbolsOut, numel(symbolsOut)/length(goldSeq), 15) * goldSeq(:, iSignal);
% end
flag = 1;
end
