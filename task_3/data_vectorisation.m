function [symbolsMatrix] = data_vectorisation(symbolsOut, nAnts, nExt, nChips)
symbolsMatrix = zeros(nAnts * nExt, length(symbolsOut) / nChips);
for iAnt = 1: nAnts
    temp = symbolsOut(:, iAnt);
    symbolsMatrix((iAnt - 1) * nExt + 1: iAnt * nExt, 1: 2: end - 1) = reshape(temp(1: length(temp) - nChips), nExt, (length(temp) - nChips) / nExt);
    symbolsMatrix((iAnt - 1) * nExt + 1: iAnt * nExt, 2: 2: end) = reshape(temp(nChips + 1: length(temp)), nExt, (length(temp) - nChips) / nExt);
end
end

