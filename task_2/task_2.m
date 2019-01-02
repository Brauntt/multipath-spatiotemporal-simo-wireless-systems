%% Initialisation
clear;close all;
% name indexes
surIndex = 3;
foreIndex = 19;
% convention
zPixel = 3;
bitInt = 8;
% maximum size of picture
widthMax = 160;
heightMax = 112;
bitsMax = widthMax * heightMax * zPixel * bitInt;
% number of signals
nSignals = 3;
% desired signal index
desiredIndex = 1;
% coefficients of the primitive polynomials
coeffs = [1 0 0 1 1; 1 1 0 0 1];
% phase shift of QPSK
phi = (surIndex + 2 * foreIndex) * pi / 180;
% directions of signals
directions = [30 0; 45 0; 20 0; 80 0; 150 0];
% number of paths per signal
nPaths = [3; 1; 1];
% path fading coefficients
fadingCoefs = [0.8; 0.4 * exp(1i * -40 / 180 * pi); 0.8 * exp(1i * 80 / 180 * pi); 0.5; 0.2];
% path delays
delays = [mod(surIndex + foreIndex, 4); 4 + mod(surIndex + foreIndex, 5); 9 + mod(surIndex + foreIndex, 6); 8; 13];
% receiver antenna array positions
array = [0 0 0];
% signal-to-noise ratio at the receiver end
snrDb = [0 40];
snr = 10 .^ (snrDb / 10);
nSnr = length(snrDb);
% minimum shift for gold sequences
shiftMin = ceil(1 + mod(surIndex + foreIndex, 12));
% maximum possible relative delay
nDelay = 2 ^ (length(coeffs) - 1) - 1;
% declaration
xPixel = zeros(nSignals, 1); yPixel = zeros(nSignals, 1); imageBits = zeros(nSignals, 1);
bitsIn = zeros(bitsMax, nSignals);
symbolsIn = zeros((2 ^ (length(coeffs) - 1) - 1) * bitsMax / 2, nSignals);
ber = zeros(nSnr, 1);
goldSeq = zeros(2 ^ (length(coeffs) - 1) - 1, nSignals);
disp(['Delays = ' num2str(delays')]);
%% Balanced gold sequence mining
% generate M-sequences
[mSeq1] = fMSeqGen(coeffs(1, :));
[mSeq2] = fMSeqGen(coeffs(2, :));
% calculate the minimum shift for balanced gold sequence
[shift] = miner(mSeq1, mSeq2, shiftMin);
%% Transmitter, channel and receiver
for iSignal = 1: nSignals
    % generate gold sequences
    goldSeq(:, iSignal) = fGoldSeq(mSeq1, mSeq2, shift + iSignal - 1);
    % declare file names
    fileName = ['pic_', num2str(iSignal), '.png'];
    % obtain the bit stream into the modulator
    [bitsIn(:, iSignal), xPixel(iSignal), yPixel(iSignal)] = fImageSource(fileName, bitsMax);
    % calculate the image size in bits
    imageBits(iSignal) = xPixel(iSignal) * yPixel(iSignal) * zPixel * bitInt;
    % modulate the signal and encode with gold sequence
    symbolsIn(:, iSignal) = fDSQPSKModulator(bitsIn(:, iSignal), goldSeq(:, iSignal), phi);
end
% show the original picture
fImageSink(bitsIn, imageBits, xPixel, yPixel);
for iSnr = 1: nSnr
    % model the channel effects in the system
    [symbolsOut] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, directions, snr(iSnr), array, nDelay);
    % estimate the delay of paths of signals
    [delayEst] = fChannelEstimation(symbolsOut{desiredIndex}, goldSeq, nPaths);
    % demodulate the received patterns
    [bitsOut] = fDSQPSKDemodulator(symbolsOut{desiredIndex}, goldSeq, phi, delayEst, nPaths, fadingCoefs);
    % display the recovered pictures
    fImageSink(bitsOut, imageBits, xPixel, yPixel, snrDb(iSnr));
    % calculate bit error rate of the desired signal
    ber(iSnr) = sum(xor(bitsOut(:, desiredIndex), bitsIn(:, desiredIndex))) / length(bitsOut);
    disp(['----------   SNR = ' num2str(snrDb(iSnr)) ' dB ----------']);
    disp(['Estimated delays = ' num2str(delayEst')]);
    disp(['Bit error rate (Source ' num2str(desiredIndex) ') = ' num2str(ber(iSnr))]);
end
% rearrange the positions of the figures
tilefigs([0 0.5 0.8 0.5]);
