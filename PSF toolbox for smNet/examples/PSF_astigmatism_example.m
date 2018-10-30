%------ Demo code for simulation of PSFs from interpolation of recoreded PSF data------------
% software requirement: Matlab R2015a or later
%                       Dipimage toolbox 2.7 or later
% system requirement:   CPU Intel Core i7
%                       32 GB RAM 
% Data format: data must be .mat file 
%
% Note:this example code demonstrates the usage of PSF_zernike and CalCRLB class, 
%      type 'help PSF_zernike' and 'help CalCRLB' in the command window for more information
%
% (C) Copyright 2018               
%     All rights reserved          
%                                 
% Author: Sheng Liu, May 2018

%% load the selected OptimPR_Ast object
basepath = pwd;
load([basepath,'/test data/astigmatimPSF_PR_result_optimAst.mat'])
probj = obj.PRobj;
%% generate PSF
psfobj = PSF_zernike(probj.PRstruct);       % create object from PSF_zernike class
Num = 17;                                   % number of PSFs
zpos = linspace(-0.8,0.8,Num)';             % z positions of the PSFs , unit is micron
psfobj.Xpos = zeros(Num,1);                 % x positions of the PSFs, unit is pixel
psfobj.Ypos = zeros(Num,1);                 % y positions of the PSFs, unit is pixel
psfobj.ZposMed = zpos;                      % z positions of the PSFs, they are the positions in the sample medium
psfobj.Boxsize = 32;                        % output size of the PSFs
psfobj.PSFsize = 128;                       % image size used for PSF generation
psfobj.Pixelsize = 0.113;                   % pixel size on the sample plane, unit is micron
psfobj.nMed = 1.35;                         % refractive index of the sample medium
psfobj.Zpos = 1;                            % stage position, unit is micron. It is a scalar when considering index mismatch aberration 

psfobj.precomputeParam();                   % generate parameters for Fourier space operation
psfobj.genPupil();                          % generate pupil function 
psfobj.genIMMPSF();                         % generate PSFs considering index mismatch aberration
psfobj.scalePSF('IMM');                     % generate OTF rescaled PSF, 'IMM': PSFs with index mismatch aberration, 'normal': PSFs without index mismatch aberration
psf = psfobj.ScaledPSFs;                    % simulated PSFs
psfobj.genfigs('IMM');                      % show simulated PSFs, 'IMM': PSFs with index mismatch aberration, 'normal': PSFs without index mismatch aberration
 
%% calculate CRLB from above PSF model
photon = 2000;                              % total photon count of a single emitter
bg = 10;                                    % background photon count of each pixel
crobj = CalCRLB(psfobj.PRstruct,'IMM');     % create object for CalCRLB class
                                            % input parameters: PRstruct,PSFtype
                                            % PSFtype - 'zernike' : PSFs without index mismatch aberration and generated from PSF_zernike class 
                                            %           'IMM'     : PSFs with index mismatch aberration and generated from PSF_zernike class 
                                            %           'DH'      : double helix PSFs   
crobj.Photon = photon.*ones(Num,1);         %
crobj.Bg = zeros(Num,1)+bg;                 %
crobj.Pixelsize = psfobj.Pixelsize;         % copy parameters from psfobj to crobj
crobj.Xpos = psfobj.Xpos;                   %
crobj.Ypos = psfobj.Ypos;                   %
crobj.Zpos = psfobj.ZposMed;                %
crobj.Boxsize = psfobj.Boxsize;             %
crobj.PSFobj.Zpos = psfobj.Zpos;            %
crobj.PSFobj.nMed = psfobj.nMed;            %
crobj.PSFobj.PSFsize = psfobj.PSFsize;      %
crobj.Deltax = 0.1;                         % increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
crobj.Deltaz = 0.01;                        % increment in z direction forcalculating first and second derivative of the objective function, unit is micron 

crobj.prepInputparam();                     % generate parameters of PSFs used for CRLB calculation
crobj.calcrlb();                            % calculate CRLB of simulated emitters, given a PSF model
crobj.genfigs();                            % generate plots of theoretical localization precision in x, y and z at z positions defined by 'crobj.Zpos'

% output: crobj.X_STD       - precision in x estimation
%         crobj.Y_STD       - precision in y estimation
%         crobj.Z_STD       - precision in z estimation
