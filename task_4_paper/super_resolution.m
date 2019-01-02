function [weight] = super_resolution(array, dirTarget, doa)
% Function: 
%   - obtain the directional gains of the superresolution beamformer
%
% InputArg(s):
%   - array: coordinates of the receiving sensors
%   - dirTarget: directions of the desired user
%   - doa: direction of arrival
%
% OutputArg(s):
%   - weight: the weight on receiving antennas
%
% Comments:
%   - superresolution beamformers maximise the SIR, which is optimum to
%   suppress interference
%
% Author & Date: Yang (i@snowztail.com) - 27 Nov 18

% obtain directions of interferences
dirInterf = setdiff(doa, dirTarget, 'rows');
% hence the manifold vector of interferences
spvInterf = spv(array, dirInterf);
% the manifold vector of the desired signal
spvTarget = spv(array, dirTarget);
% weights of super-resolution beamformer
weight = fpoc(spvInterf) * spvTarget;
end

