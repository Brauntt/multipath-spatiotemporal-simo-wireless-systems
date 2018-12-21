function [bitsIn, xPixel, yPixel] = fImageSource(fileName, bitsMax)
% Function:
%   - reads an image file with AxB pixels and produces a column vector of
%  bits of length Q=AxBx3x8 where 3 represents the R, G and B matrices used
%  to represent the image and 8 represents an 8 bit integer. If bitsMax>Q
%  then the vector is padded at the bottom with zeros.
%
% InputArg(s):
%   - fileName: the file name of the image
%   - bitsMax: the upper bound of the picture size in bits
%
% OutputArg(s):
%   - bitsIn: bitsMax bits representing the image
%   - xPixel: number of pixels in image in x dimension
%   - yPixel: number of pixels in image in y dimension
%
% Comments:
%   - the input file should have a smaller resolution than A x B
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% bit integer
bitInt = 8;
% obtain the RGB of the picture
rgb = imread(fileName);
% obtain the pixels in x and y dimension
[xPixel, yPixel, ~] = size(rgb);
% convert RGB to binary
rgbBin = de2bi(rgb, bitInt);
% reshape as a stream
bitsIn = reshape(rgbBin, numel(rgbBin), 1);
% pad zeros at the end
if length(bitsIn) < bitsMax
    bitsIn(bitsMax) = 0;
end
% convert the variable type to double
bitsIn = double(bitsIn);
end
