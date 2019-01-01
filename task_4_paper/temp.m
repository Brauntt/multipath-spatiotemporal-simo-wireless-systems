function [objTemporalSmooth] = temp(nSubVects, lenSubVect, nChips, obj)

nSegs = length(obj) / (2 * nChips);
for iSegRow = 1: nSegs
    for iSegCol = 1: nSegs
       objSeg = obj((iSegRow - 1) * 2 * nChips + 1: iSegRow * 2 * nChips, (iSegCol - 1) * 2 * nChips + 1: iSegCol * 2 * nChips);
       objSub = 0;
       for iSubVect = 1: nSubVects
           objSub = objSub + objSeg(iSubVect: iSubVect + lenSubVect - 1, iSubVect: iSubVect + lenSubVect - 1);
       end
       objSub = objSub / nSubVects;
       objTemporalSmooth((iSegRow - 1) * lenSubVect + 1: iSegRow * lenSubVect, (iSegCol - 1) * lenSubVect + 1: iSegCol * lenSubVect) = objSub;
    end
end
end
