% NAME, GROUP (EE4/MSc), 2010, Imperial College.
% DATE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display the received image by converting bits back into R, B and G
% matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% bitsIn (Px1 Integers) = P demodulated bits of 1's and 0's
% Q (Integer) = Number of bits in the image
% x (Integer) = Number of pixels in image in x dimension
% y (Integer) = Number of pixels in image in y dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fImageSink(bitsOut, Q, x, y, snrDb)
nSignals = length(Q);
rgb = cell(nSignals, 1);
figure;
for iSignal = 1: nSignals
    bits = bitsOut(1: Q(iSignal), iSignal);
    bits = uint8(bits);
    rgbBin = reshape(bits, length(bits)/8, 8);
    rgb{iSignal} = reshape(bi2de(rgbBin), x(iSignal), y(iSignal), 3);
    subplot(nSignals, 1, iSignal);
    imshow(rgb{iSignal});
    title(['Image of user ', num2str(iSignal), ' (snr = ', num2str(snrDb) ' dB)']);
end
end
