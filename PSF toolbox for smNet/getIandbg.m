% GETIANDBG    estimates total photon and background from single molecule subregions
%
% SYNOPSIS:
%   [I,bg] = getIandbg(data)
%
% INPUTS:
%   data
%       a stack of subregions with a single molecule near the center of each subregion
%
% OUTPUTS:
%   I
%       estimations of total photon counts of each subregion
%   bg
%       estimations of background photon counts of each subregion
%
% (C) Copyright 2017                Huang Lab, Weldon School of Biomedical Engineering, 
%     All rights reserved           Purdue University, West Lafayette, IN, USA
%
%                                   
% Author: Sheng Liu, July 2018


function [I,bg] = getIandbg(data)
sz = size(data);
if numel(sz)<3
    sz(3) = 1;
end
rect = ones(sz(1),sz(2));
rect(2:end-1,2:end-1) = 0;
mask = rect ~= 0;
dataxz = reshape(data,sz(1)*sz(1),sz(3));
bg = ones(1,1,sz(3));
bg(1,1,:) = median(dataxz(mask(:),:),1)';
BG = repmat(bg.*0.9,[sz(1),sz(2),1]);
I = squeeze(sum(sum(data-BG,1),2));
I(I<200) = 200;
bg = squeeze(bg);
end