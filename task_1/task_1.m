clear;close;
surIndex = 26;
foreIndex = 25;
widthMax = 160;
heightMax = 112;
nSignals = 3;
coeffs = [1 0 0 1 1; 1 1 0 0 1];
phi = (surIndex + 2 * foreIndex) * pi / 180;
% p = 1e6;
p = widthMax * heightMax * 3 * 8;
doa = [30 0; 90 0; 150 0];
array = [0 0 0];
snrDb = [0 40];
snr = 10 .^ (snrDb / 10);
x = zeros(nSignals, 1); y = zeros(nSignals, 1); Q = zeros(nSignals, 1);
%% Gold sequence
symbolsIn = zeros((2 ^ (length(coeffs) - 1) - 1) * p / 2, nSignals);
nPaths = [1; 1; 1];
delays = [5; 7; 12];
fadingCoefs = [0.4; 0.7; 0.2];
[mSeq1] = fMSeqGen(coeffs(1, :));
[mSeq2] = fMSeqGen(coeffs(2, :));
shiftMin = ceil(1 + mod(surIndex + foreIndex, 12));
[shift] = miner(mSeq1, mSeq2, shiftMin);
goldSeq = zeros(2 ^ (length(coeffs) - 1) - 1, nSignals);
%% Signal generation
nSnr = length(snrDb);
for iSnr = 1: nSnr
    for iSignal = 1: nSignals
        goldSeq(:, iSignal) = fGoldSeq(mSeq1, mSeq2, shift + iSignal - 1);
        fileName = ['pic_', num2str(iSignal), '.png'];
        % fileName = [num2str(iSignal), '.jpg'];
        [bitsIn, x(iSignal), y(iSignal)] = fImageSource(fileName, p);
        Q(iSignal) = x(iSignal) * y(iSignal) * 3 * 8;
        % fImageSink(bitsIn, Q, x, y);
        symbolsIn(:, iSignal) = fDSQPSKModulator(bitsIn, goldSeq(:, iSignal), phi);
    end
    [symbolsOut] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, doa, snr(iSnr), array, goldSeq);
    [bitsOut] = fDSQPSKDemodulator(symbolsOut, goldSeq, phi);
    fImageSink(bitsOut, Q, x, y, snrDb(iSnr));
end
flag = 1;
