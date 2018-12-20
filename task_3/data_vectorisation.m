function [symbolsMatrix] = data_vectorisation(symbolsOut, nAnts, nChips)
if mod(length(symbolsOut), nChips)
    symbolsOut(length(symbolsOut) + mod(length(symbolsOut), nChips), end) = 0;
end
symbolsMatrix = zeros(nAnts * 2 * nChips, length(symbolsOut) / nChips);
for iAnt = 1: nAnts
    temp = symbolsOut(:, iAnt);
    symbolsMatrix((iAnt - 1) * 2 * nChips + 1: iAnt * 2 * nChips, 1: 2: end - 1) = reshape(temp(1: length(temp) - nChips), 2 * nChips, (length(temp) - nChips) / (2 * nChips));
    symbolsMatrix((iAnt - 1) * 2 * nChips + 1: iAnt * 2 * nChips, 2: 2: end) = reshape(temp(nChips + 1: length(temp)), 2 * nChips, (length(temp) - nChips) / (2 * nChips));
end
end
