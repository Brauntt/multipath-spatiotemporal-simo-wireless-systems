function [symbolsOut, desiredNoisePower] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, directions, snr, array, nDelay, desiredIndex)
% Function:
%   - Model the channel effects in the system
%
% InputArg(s):
%   - nPaths: number of paths for each source in the system. For example,
%  if 3 sources with 1, 3 and 2 paths respectively then nPaths = [1;3;2]
%   - symbolsIn: input symbols of the channel
%   - delays: delays for each path in the system starting with source 1
%   - fadingCoefs: fading coefficients for each path in the system
%   - directions: direction of arrival for each source in the system in the
%  form [Azimuth, Elevation]
%   - snr: signal-to-noise ratio at the receiver end
%   - array: array locations in half unit wavelength. If no array then
%  should be [0,0,0]
%   - nDelay: maximum possible relative delay
%   - desiredIndex: desired signal index
%
% OutputArg(s):
%   - symbolsOut: channel symbol chips received from each antenna
%   - desiredNoisePower: noise power added to desired signal
%
% Comments:
%   - each antenna receive different copies of signals
%   - the antennas regard a particular signal as desired and the noise
%  should be calculated with the target signal power and the SNR. The
%  noise should be added in channel but it is yet not clear which signal
%  is desired. Therefore, noise is calculated for each signal and added to
%  the ideal output individually, which produces multiple output symbol
%  streams with noise level corresponding to different target signals. The
%  receiver end only need to specify which one is the desired signal and
%  use the corresponding output stream.
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% number of signals
nSignals = length(nPaths);
% expand the symbol streams by maximum possible relative delay
symbolsIn(length(symbolsIn) + nDelay, nSignals) = 0;
% path counters
pathCounter = 1; pathCounterPrev = 1;
% symbol streams of signals
symbolsSignal = zeros(size(symbolsIn));
% symbol streams of paths
symbolsPath = zeros(length(symbolsIn), sum(nPaths));
% noise should be calculated individually different target signals
noise = cell(nSignals, 1);
% the output streams for different target signals, determined by varied
% noise levels
symbolsOut = cell(nSignals, 1);
% gain of elements on users directions
spvSources = spv(array, directions);
for iSignal = 1: nSignals
    % symbol streams of paths of the current signal
    symbolsSub = zeros(length(symbolsIn), nPaths(iSignal));
    for iPath = 1: nPaths(iSignal)
        % obtain the current path stream by delay and scaling
        symbolsSub(:, iPath) = fadingCoefs(pathCounter) * circshift(symbolsIn(:, iSignal), delays(pathCounter));
        % store the current path stream
        symbolsPath(:, pathCounter) = symbolsSub(:, iPath);
        % update the path counter
        pathCounter = pathCounter + 1;
    end
    % sum the path streams to get signal stream
    symbolsSignal(:, iSignal) = sum(symbolsSub, 2);
    % assume the current signal as desired
    symbolsDesired = (spv(array, directions(pathCounterPrev: pathCounter - 1, :)) * symbolsSub.').';
    % update the previous path counter
    pathCounterPrev = pathCounter;
    % calculate the desired signal power
    powerSignal = mean(sum(abs(symbolsDesired).^2) / length(symbolsDesired));
    % hence the noise power
    powerNoise = powerSignal / snr;
    % obtain the noise power of the desired signal
    if iSignal == desiredIndex
       desiredNoisePower = powerNoise;
    end
    % then generate noise accordingly
    noise{iSignal} =  sqrt(powerNoise / 2) * (randn(size(symbolsDesired)) + 1i * randn(size(symbolsDesired)));
end
for iSignal = 1: nSignals
    % produces multiple output symbol streams with noise level based on
    % different target signals
    symbolsOut{iSignal} = (spvSources * symbolsPath.').' + noise{iSignal};
end
end
