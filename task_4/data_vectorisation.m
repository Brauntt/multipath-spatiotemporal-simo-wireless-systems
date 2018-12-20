function [symbolsMatrix] = data_vectorisation(signalSample, nAnts, nChips)
if mod(length(signalSample), nChips)
    signalSample(length(signalSample) + mod(length(signalSample), nChips), end) = 0;
end
symbolsMatrix = zeros(nAnts * 2 * nChips, length(signalSample) / nChips);
for iAnt = 1: nAnts
    temp = signalSample(:, iAnt);
    symbolsMatrix((iAnt - 1) * 2 * nChips + 1: iAnt * 2 * nChips, 1: 2: end - 1) = reshape(temp(1: length(temp) - nChips), 2 * nChips, (length(temp) - nChips) / (2 * nChips));
    symbolsMatrix((iAnt - 1) * 2 * nChips + 1: iAnt * 2 * nChips, 2: 2: end) = reshape(temp(nChips + 1: length(temp)), 2 * nChips, (length(temp) - nChips) / (2 * nChips));
end
end
