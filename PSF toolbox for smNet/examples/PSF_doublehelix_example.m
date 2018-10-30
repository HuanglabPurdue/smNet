%------ Demo code for simulation of double-helix PSFs------------
% software requirement: Matlab R2015a or later
%                       Dipimage toolbox 2.7 or later
% system requirement:   CPU Intel Core i7
%                       32 GB RAM 
% Note:this example code demonstrates the usage of PSF_DH and CalCRLB class, 
%      type 'help PSF_DH' and 'help CalCRLB' in the command window for more information

% (C) Copyright 2018               
%     All rights reserved           
%
% Author: Sheng Liu, May 2018

%% create PRstruct and setup parameters
PRstruct.NA = 1.49;                                 % numerical aperture of the objective lens
PRstruct.Lambda = 0.69;                             % center wavelength of the emission band pass filter, unit is micron
PRstruct.RefractiveIndex = 1.52;                    % refractive index of the immersion medium
%% generate PSF
psfobj = PSF_DH(PRstruct);                          % create object from PSF_DH class
Num = 21;                                           % number of PSFs
zpos = linspace(-1,1,Num)';                         % z positions of the PSFs , unit is micron
psfobj.Xpos = zeros(Num,1);                         % x positions of the PSFs, unit is pixel
psfobj.Ypos = zeros(Num,1);                         % y positions of the PSFs, unit is pixel
psfobj.Zpos = zpos;                                 % 
psfobj.Boxsize = 32;                                % output size of the PSFs
psfobj.Pixelsize = 0.1;                             % pixel size on the sample plane, unit is micron
psfobj.Magnification = 100;                         % magnification of the optical system
psfobj.PSFsize = 512;                               % image size used for PSF generation
psfobj.Zlim = [-1.2,1.2];                           % axial range for optimization of the pupil phase

psfobj.precomputeParam();                           % generate parameters for Fourier space operation
psfobj.genPupil();                                  % generate high efficiency pupil function, takes ~1 minute
psfobj.genPSF();                                    % generate PSFs
psf = psfobj.ScaledPSFs;                            % simulated PSF, they are without OTF rescale, but here used the same name for consistency
psfobj.genfigs();                                   % show simulated PSFs

%% calculate CRLB from above PSF model
photon = 2000;                                      % total photon count of a single emitter                            
bg = 5;                                             % background photon count of each pixel
crobj = CalCRLB(PRstruct,'DH');                     % create object for CalCRLB class
                                                    % input parameters: PRstruct,PSFtype
                                                    % PSFtype - 'zernike' : PSFs without index mismatch aberration and generated from PSF_zernike class 
                                                    %           'IMM'     : PSFs with index mismatch aberration and generated from PSF_zernike class 
                                                    %           'DH'      : double helix PSFs   
crobj.Photon = photon.*ones(Num,1);                 %
crobj.Bg = zeros(Num,1)+bg;                         %
crobj.Pixelsize = psfobj.Pixelsize;                 % --------copy parameters from psfobj to crobj---------
crobj.Xpos = psfobj.Xpos;                           %
crobj.Ypos = psfobj.Ypos;                           %
crobj.Zpos = psfobj.Zpos;                           %
crobj.Boxsize = psfobj.Boxsize;                     %
crobj.PSFobj.Magnification = psfobj.Magnification;  %
crobj.PSFobj.PSFsize = psfobj.PSFsize;              %
crobj.PSFobj.Zlim = psfobj.Zlim;                    % -----------------------------------------------------
crobj.Deltax = 0.1;                                 % increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
crobj.Deltaz = 0.01;                                % increment in z direction forcalculating first and second derivative of the objective function, unit is micron 
 
crobj.prepInputparam();                             % generate parameters of PSFs used for CRLB calculation
crobj.calcrlb();                                    % calculate CRLB of simulated emitters, given a PSF model
crobj.genfigs();                                    % generate plots of theoretical localization precision in x, y and z at z positions defined by 'crobj.Zpos'

% output: crobj.X_STD       - precision in x estimation
%         crobj.Y_STD       - precision in y estimation
%         crobj.Z_STD       - precision in z estimation

