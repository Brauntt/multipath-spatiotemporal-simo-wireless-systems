function display_text(bitsOut, nChars)
% Function:
%   - convert the demodulated bits to text
%
% InputArg(s):
%   - bitsOut: demodulated bit stream
%   - nChars: number of characters of the text
%
% OutputArg(s):
%   - none (display the text in the command window)
%
% Comments:
%   - text size should be known
%
% Author & Date: Yang (i@snowztail.com) - 1 Jan 19

% convert the bits to desired class and reshape
bitsText = uint8(reshape(bitsOut, length(bitsOut) / nChars, nChars));
% then to decimal
bitsText = bi2de(bitsText.', 'left-msb')';
% display the text
disp(char(bitsText));
end

