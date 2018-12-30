function [symbolsMatrix] = data_vectorisation(symbolsOut, nAnts, nChips)
% Function:
%   - perform data vectorisation for symbol stream
%
% InputArg(s):
%   - symbolsOut: channel symbol chips received
%   - nAnts: number of antennas
%   - nChips: chip length
%
% OutputArg(s):
%   - symbolsMatrix: 2-D symbol matrix
%
% Comments:
%   - realise by assigning odd and even columns
%
% Author & Date: Yang (i@snowztail.com) - 21 Dec 18

% ensure the stream length is integer times of chip length by padding zeros
if mod(length(symbolsOut), nChips)
    symbolsOut(length(symbolsOut) + mod(length(symbolsOut), nChips), end) = 0;
end
% declare size of dymbol matrix
symbolsMatrix = zeros(nAnts * 2 * nChips, (length(symbolsOut) - nChips) / nChips);
for iAnt = 1: nAnts
    % symbols on the current antenna
    symbolsAnt = symbolsOut(:, iAnt);
    % odd columns of symbol matrix
    symbolsMatrix((iAnt - 1) * 2 * nChips + 1: iAnt * 2 * nChips, 1: 2: end - 1) = reshape(symbolsAnt(1: length(symbolsAnt) - nChips), 2 * nChips, (length(symbolsAnt) - nChips) / (2 * nChips));
    % even columns of symbol matrix
    symbolsMatrix((iAnt - 1) * 2 * nChips + 1: iAnt * 2 * nChips, 2: 2: end) = reshape(symbolsAnt(nChips + 1: length(symbolsAnt)), 2 * nChips, (length(symbolsAnt) - nChips) / (2 * nChips));
end
end
