function [objSmooth] = spatial_smoothing(lenSubVect, nAnts, obj)
objSmooth = 0;
lenSmooth = lenSubVect * nAnts;
nSmoothes = length(obj) - lenSmooth + 1;
for iSmooth = 1: nSmoothes
    objSmooth = objSmooth + obj(iSmooth: iSmooth + lenSmooth - 1, iSmooth: iSmooth + lenSmooth - 1);
end
objSmooth = objSmooth / nSmoothes;
end

