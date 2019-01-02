%% Initialisation
clear; close all;
% load data
load('yz16718.mat');
% coefficients of the primitive polynomials
coeffs = [1 0 0 1 0 1; 1 0 1 1 1 1];
% shift of the second M-sequence to create gold sequence
shift = phase_shift;
% phase shift of QPSK
phi = phi_mod * pi / 180;
% path fading coefficients
fadingCoefs = Beta_1;
% number of characters of the text
nChars = 60;
% number of paths of the desired signal
nPaths = length(fadingCoefs);
% chip length
nChips = 2 ^ (length(coeffs) - 1) - 1;
% signal samples
signalSample = Xmatrix.';
% number of snapshots and antennas
[nSnapshots, nAnts] = size(signalSample);
% initial phase of antenna positions
initPhase = 30 / 180 * pi;
% receiver antenna array positions normalised to half wavelengths
array = zeros(nAnts, 3);
for iAnt = 1: nAnts
    % array positions
    array(iAnt, :) = [cos(initPhase + (iAnt - 1) * 2 * pi / 5), sin(initPhase + (iAnt - 1) * 2 * pi / 5), 0];
end
%% Gold sequence generation
% generate M-sequences
[mSeq1] = fMSeqGen(coeffs(1, :));
[mSeq2] = fMSeqGen(coeffs(2, :));
% gold sequence
goldSeq = fGoldSeq(mSeq1, mSeq2, shift);
%% Detection and estimation
% data vectorisation
[symbolsMatrix] = data_vectorisation(signalSample, nAnts, nChips);
% estimate the delay and DOA of paths of signals
[doaEst, delayEst] = fChannelEstimation(array, signalSample, goldSeq, nPaths);
%% Demodulation
% obtain the weights of the spatiotemporal rake beamformer
[weightStRake] = spatiotemporal_rake(array, doaEst, delayEst, goldSeq, nPaths, fadingCoefs);
% demodulate the signal
[bitsOut] = fDSQPSKDemodulator(symbolsMatrix, weightStRake, goldSeq, phi);
% display the text
display_text(bitsOut, nChars);
disp(['Estimated delays = ' num2str(delayEst')]);
disp(['Estimated DOAs = ' num2str(reshape(doaEst', 1, numel(doaEst)))]);
% rearrange the positions of the figures
tilefigs([0 0.5 0.8 0.5]);
