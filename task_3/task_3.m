clear; close;
nAnts = 5;
initPhase = 30 / 180 * pi;
array = zeros(nAnts, 3);
for iAnt = 1: nAnts
    array(iAnt, :) = [cos(initPhase + (iAnt - 1) * 2 * pi / 5), sin(initPhase + (iAnt - 1) * 2 * pi / 5), 0];
end
azimuth = 0: 180;
elevation = 0;
directions = [30 0; 45 0; 20 0; 80 0; 150 0];
snrDb = [40 40];
nSnr = length(snrDb);
surIndex = 26;
foreIndex = 25;
widthMax = 160;
heightMax = 112;
nSignals = 3;
coeffs = [1 0 0 1 1; 1 1 0 0 1];
phi = (surIndex + 2 * foreIndex) * pi / 180;
% p = 1e6;
p = widthMax * heightMax * 3 * 8;
varNoise = 10 .^ (-snrDb / 10);
snr = 10 .^ (snrDb / 10);
x = zeros(nSignals, 1); y = zeros(nSignals, 1); Q = zeros(nSignals, 1);
bitsIn = zeros(p, nSignals);
%% Gold sequence
symbolsIn = zeros((2 ^ (length(coeffs) - 1) - 1) * p / 2, nSignals);
nPaths = [3; 1; 1];
delays = [mod(surIndex + foreIndex, 4); 4 + mod(surIndex + foreIndex, 5); 9 + mod(surIndex + foreIndex, 6); 8; 13];
fadingCoefs = [0.8; 0.4 * exp(1i * -40 / 180 * pi); 0.8 * exp(1i * 80 / 180 * pi); 0.5; 0.2];
[mSeq1] = fMSeqGen(coeffs(1, :));
[mSeq2] = fMSeqGen(coeffs(2, :));
shiftMin = ceil(1 + mod(surIndex + foreIndex, 12));
[shift] = miner(mSeq1, mSeq2, shiftMin);
goldSeq = zeros(2 ^ (length(coeffs) - 1) - 1, nSignals);
nChips = length(goldSeq);
nExt = 2 * nChips;
% shiftMatrix = [zeros(1, nExt); eye(nExt - 1) zeros(nExt - 1, 1)];
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
% fImageSink(bitsIn, Q, x, y);
symbolsMatrix = zeros(nAnts * nExt, length(symbolsIn) / nChips);
for iSnr = 1: nSnr
    [symbolsOut] = fChannel(nPaths, symbolsIn, delays, fadingCoefs, directions, snr(iSnr), array, goldSeq);
    for iAnt = 1: nAnts
        temp = symbolsOut(:, iAnt);
        symbolsMatrix((iAnt - 1) * nExt + 1: iAnt * nExt, 1: 2: end - 1) = reshape(temp(1: length(temp) - nChips), nExt, (length(temp) - nChips) / nExt);
        symbolsMatrix((iAnt - 1) * nExt + 1: iAnt * nExt, 2: 2: end) = reshape(temp(nChips + 1: length(temp)), nExt, (length(temp) - nChips) / nExt);
    end
    covSymbol = symbolsMatrix * symbolsMatrix' / length(symbolsMatrix);
    doa = music(array, covSymbol, goldSeq);
end
