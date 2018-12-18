clear;close all;
surIndex = 26;
foreIndex = 25;
widthMax = 160;
heightMax = 112;
nSignals = 3;
coeffs = [1 0 0 1 1; 1 1 0 0 1];
directions = [30 0; 45 0; 20 0; 80 0; 150 0];
phi = (surIndex + 2 * foreIndex) * pi / 180;
% p = 1e6;
p = widthMax * heightMax * 3 * 8;
array = [0 0 0];
snrDb = [10 40];
snr = 10 .^ (snrDb / 10);
nPaths = [3; 1; 1];
nSnr = length(snrDb);
bitsIn = zeros(p, nSignals);
% ber = zeros(nSnr, nSignals);
ber = zeros(nSnr, 1);
delays = [mod(surIndex + foreIndex, 4); 4 + mod(surIndex + foreIndex, 5); 9 + mod(surIndex + foreIndex, 6); 8; 13]
fadingCoefs = [0.8; 0.4 * exp(1i * -40 / 180 * pi); 0.8 * exp(1i * 80 / 180 * pi); 0.5; 0.2];
% varNoise = 10 .^ (-snrDb / 10);
x = zeros(nSignals, 1); y = zeros(nSignals, 1); Q = zeros(nSignals, 1);
desiredUserIndex = 1;
%% Gold sequence
symbolsIn = zeros((2 ^ (length(coeffs) - 1) - 1) * p / 2, nSignals);
[mSeq1] = fMSeqGen(coeffs(1, :));
[mSeq2] = fMSeqGen(coeffs(2, :));
shiftMin = ceil(1 + mod(surIndex + foreIndex, 12));
[shift] = miner(mSeq1, mSeq2, shiftMin);
% shift = 6;
goldSeq = zeros(2 ^ (length(coeffs) - 1) - 1, nSignals);
%% Signal generation
for iSignal = 1: nSignals
    goldSeq(:, iSignal) = fGoldSeq(mSeq1, mSeq2, shift + iSignal - 1);
%     fileName = ['t', num2str(iSignal), '.jpg'];
    fileName = ['pic_', num2str(iSignal), '.png'];
%     fileName = [num2str(iSignal), '.jpg'];
    [bitsIn(:, iSignal), x(iSignal), y(iSignal)] = fImageSource(fileName, p);
    Q(iSignal) = x(iSignal) * y(iSignal) * 3 * 8;
    % fImageSink(bitsIn, Q, x, y);
    symbolsIn(:, iSignal) = fDSQPSKModulator(bitsIn(:, iSignal), goldSeq(:, iSignal), phi);
end
for iSnr = 1: nSnr
    [symbolsOut] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, directions, snr(iSnr), array, goldSeq);
    % desired user index is 1
    [delayEst] = fChannelEstimation(symbolsOut{desiredUserIndex}, goldSeq, nPaths)
%     delayEst = delays;
    [bitsOut] = fDSQPSKDemodulator(symbolsOut{desiredUserIndex}, goldSeq, phi, delayEst, nPaths, fadingCoefs);
    fImageSink(bitsOut, Q, x, y, snrDb(iSnr));
%     ber(iSnr, :) = sum(xor(bitsOut, bitsIn));
    ber(iSnr) = sum(xor(bitsOut(:, 1), bitsIn(:, 1))) / length(bitsOut)
end
fImageSink(bitsIn, Q, x, y);
