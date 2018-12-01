% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Models the channel effects in the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% paths (Mx1 Integers) = Number of paths for each source in the system.
% For example, if 3 sources with 1, 3 and 2 paths respectively then
% paths=[1;3;2]
% symbolsIn (MxR Complex) = Signals being transmitted in the channel
% delay (Cx1 Integers) = Delay for each path in the system starting with
% source 1
% beta (Cx1 Integers) = Fading Coefficient for each path in the system
% starting with source 1
% DOA = Direction of Arrival for each source in the system in the form
% [Azimuth, Elevation]
% SNR = Signal to Noise Ratio in dB
% array = Array locations in half unit wavelength. If no array then should
% be [0,0,0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% symbolsOut (FxN Complex) = F channel symbol chips received from each antenna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [symbolsOut] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, doa, snr, array, goldSeq)
[nDelays, nSignals] = size(goldSeq);
symbolsIn(length(symbolsIn) + nDelays, nSignals) = 0;
symbolsUser = zeros(size(symbolsIn));
pathCounter = 1;
for iSignal = 1: nSignals
    symbolsPath = zeros(length(symbolsIn), nPaths(iSignal));
    for iPath = 1: nPaths(iSignal)
        symbolsPath(:, iPath) = fadingCoefs(pathCounter) * circshift(symbolsIn(:, iSignal), delays(pathCounter));
        noise = 1 / sqrt(2 * snr) * (randn(length(symbolsIn), 1) + 1i * randn(length(symbolsIn), 1));
        symbolsPath(:, iPath) = symbolsPath(:, iPath) + noise;
        pathCounter = pathCounter + 1;
    end
    symbolsUser(:, iSignal) = sum(symbolsPath, 2);
end
symbolsOut = sum(symbolsUser, 2);
end
