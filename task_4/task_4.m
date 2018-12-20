clear; close;
load('yz16718.mat');
coeffs = [1 0 0 1 0 1; 1 0 1 1 1 1];
shift = phase_shift;
phi = phi_mod * pi / 180;
fadingCoefs = Beta_1;
nPaths = length(fadingCoefs);
signalSample = Xmatrix';
[nSnapshots, nAnts] = size(signalSample);
initPhase = 30 / 180 * pi;
array = zeros(nAnts, 3);
for iAnt = 1: nAnts
    array(iAnt, :) = [cos(initPhase + (iAnt - 1) * 2 * pi / 5), sin(initPhase + (iAnt - 1) * 2 * pi / 5), 0];
end
%% Gold sequence generation
[mSeq1] = fMSeqGen(coeffs(1, :));
[mSeq2] = fMSeqGen(coeffs(2, :));
goldSeq = fGoldSeq(mSeq1, mSeq2, shift);
nChips = length(goldSeq);
%% Detection and estimation
[symbolsMatrix] = data_vectorisation(signalSample, nAnts, nChips);
covSample = symbolsMatrix * symbolsMatrix' / length(symbolsMatrix);
[nSources] = detector_mdl(nSnapshots, covSample);
[doaEst, delayEst] = music(array, symbolsMatrix, covSample, goldSeq, nPaths);
