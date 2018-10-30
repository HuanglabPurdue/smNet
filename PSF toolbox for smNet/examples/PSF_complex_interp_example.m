%------ Demo code for simulation of PSFs from interpolation of experimental recoreded PSF data------------
% software requirement: Matlab R2015a or later
%                       Dipimage toolbox 2.7 or later
% system requirement:   CPU Intel Core i7
%                       32 GB RAM 
% Data format: data must be .mat file, containing a R by R by N matrix, R
%              is the x and y dimension of the image, N is the number of
%              axial positions, with a step size of 10~40 nm
% Note:this example code demonstrates the usage of psf_interp class, 
%      type 'help psf_interp' in the command window for more information
%
% (C) Copyright 2018               
%     All rights reserved           
%
% Author: Sheng Liu, May 2018

%% load recorded sample psf data
basepath = pwd;
F = load([basepath,'/test data/complexPSF_samplepsf.mat']);
namei = fields(F);
psfdata = F.(namei{1});
psfdata = squeeze(sum(psfdata,3));
offset = 100;                       % camera offset of recorded PSF data, unit is ADU
gain = 2;                           % camera gain of recorded PSF data
psfdata = (psfdata-offset)./gain;   % convert ADU count to photon count
%% generate PSFs
psfobj = PSF_interp(psfdata);       % create object from PSF_interp class 
Num = 25;                           % number of PSFs
zpos = linspace(-1,1,Num)';         % z positions of the PSFs , unit is micron
psfobj.Xpos = zeros(Num,1);         % x positions of the PSFs, unit is pixel
psfobj.Ypos = zeros(Num,1);         % y positions of the PSFs, unit is pixel
psfobj.Zpos = zpos;                 % 
psfobj.Boxsize = 32;                % output size of the PSFs
psfobj.Pixelsize = 0.129;           % pixel size on the sample plane, unit is micron
psfobj.PSFsize = 50;                % size of the cropped region from the recorded PSF data, this defines the x,y range of the sample psf, usually should be 1.5 times larger than psfobj.Boxsize
psfobj.Zrange = [-1.5,1.5];         % z range of sample psf data, unit is micron
 
psfobj.genSamplePSF();              % generate sample PSF from recorded PSF data, which consists of the sample points for interpolation
psfobj.genPSF();                    % generate PSFs
psf = psfobj.ScaledPSFs;            % interpolated PSF, they are without OTF rescale, but here used the same name for consistency
psfobj.genfigs();                   % show simulated PSFs

