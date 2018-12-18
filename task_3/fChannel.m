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

function [symbolsOut] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, directions, snr, array, goldSeq)
[nRelativeDelays, nSignals] = size(goldSeq);
symbolsIn(length(symbolsIn) + nRelativeDelays, nSignals) = 0;
symbolsUser = zeros(size(symbolsIn));
symbolsAll = zeros(length(symbolsIn), sum(nPaths));
% noise = zeros(size(symbolsIn));
pathCounter = 1;
startCounter = 1;
% nAnts = length(array);
noise = cell(nSignals, 1);
symbolsOut = cell(nSignals, 1);
% gain of elements on users directions
spvSources = spv(array, directions);
for iSignal = 1: nSignals
    symbolsPath = zeros(length(symbolsIn), nPaths(iSignal));
    for iPath = 1: nPaths(iSignal)
        symbolsPath(:, iPath) = fadingCoefs(pathCounter) * circshift(symbolsIn(:, iSignal), delays(pathCounter));
%         noise = 0;
%         noise = fadingCoefs(pathCounter) * sqrt(varNoise) / sqrt(2) * (randn(length(symbolsIn), 1) + 1i * randn(length(symbolsIn), 1));
%         noise = sqrt(varNoise) * (randn(length(symbolsIn), 1) + 1i * randn(length(symbolsIn), 1));
%         symbolsPath(:, iPath) = symbolsPath(:, iPath) + noise;
        symbolsAll(:, pathCounter) = symbolsPath(:, iPath);
        pathCounter = pathCounter + 1;
    end
    symbolsUser(:, iSignal) = sum(symbolsPath, 2);
    symbolsIdeal = (spv(array, directions(startCounter: pathCounter - 1, :)) * symbolsPath.').';
    startCounter = pathCounter;
    powerSignal = sum(abs(symbolsIdeal).^2) / length(symbolsIdeal);
    powerNoise = powerSignal / snr;
    noise{iSignal} = (randn(length(symbolsIn), 1) + 1i * randn(length(symbolsIn), 1)) * powerNoise;
end
% noise = sqrt(varNoise) * (randn(length(symbolsIn), 1) + 1i * randn(length(symbolsIn), 1));
for iSignal = 1: nSignals
    symbolsOut{iSignal} = (spvSources * symbolsAll.').' + noise{iSignal};
end
% symbolsOut = (spvSources * symbolsAll.').' + noise{1};
% symbolsOut = sum(symbolsUser, 2);
end
