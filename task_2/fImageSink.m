function fImageSink(bitsOut, imageBits, xPixel, yPixel, snrDb)
% Function:
%   - display the received image by converting bits back into R, G and B
%  matrices
%
% InputArg(s):
%   - bitsOut: demodulated bits
%   - imageBits: number of bits in the image
%   - xPixel: number of pixels in image in x dimension
%   - yPixel: number of pixels in image in y dimension
%   - snrDb: signal-to-noise ratio in dB
%
% OutputArg(s):
%   - None
%
% Comments:
%   - this function displays the image from demodulated RGB bits
%   - if no SNR is specified, the result is regarded as the original image
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% convention
zPixel = 3;
bitInt = 8;
% number of signals
nSignals = length(imageBits);
% RGB of each signal
rgb = cell(nSignals, 1);
% plot the picture
figure;
for iSignal = 1: nSignals
    % retrieve the meaningful bits
    bits = bitsOut(1: imageBits(iSignal), iSignal);
    % convert the stream to uint8 format
    bits = uint8(bits);
    % reshape the bits stream to standard format
    rgbBin = reshape(bits, length(bits) / bitInt, bitInt);
    rgb{iSignal} = reshape(bi2de(rgbBin), xPixel(iSignal), yPixel(iSignal), zPixel);
    % display the recovered pictures
    subplot(nSignals, 1, iSignal);
    imshow(rgb{iSignal});
    if nargin == 4
        title(['Original Image of Source ', num2str(iSignal)]);
    end
    if nargin == 5
        title(['Recovered Image of Source ', num2str(iSignal), ' (SNR = ', num2str(snrDb) ' dB)']);
    end
end
end
