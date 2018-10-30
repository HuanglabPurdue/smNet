%------ Demo code for simulation of complex PSFs------------
% software requirement: Matlab R2015a or later
%                       Dipimage toolbox 2.7 or later
% system requirement:   CPU Intel Core i7
%                       32 GB RAM 
% Note:this example code demonstrates the usage of PSF_zernike and CalCRLB class, 
%      type 'help PSF_zernike' and 'help CalCRLB' in the command window for more information

% (C) Copyright 2018               
%     All rights reserved          
%                                   
% Author: Sheng Liu, May 2018

%% create PRstruct and setup parameters
R = 128;                                % image size used for PSF generation
phaseZ = zeros(1,25);                   % setup Zernike coefficients for both the magnitude and phase part of the pupil function
magZ = zeros(1,25);                     %
phaseZ([5 16 19]) = [3 4 3];            % add aberrations to pupil phase  
magZ(1) = 1;                            % 
PRstruct.Zernike_phase = phaseZ;        % 
PRstruct.Zernike_mag = magZ;            %
PRstruct.NA = 1.4;                      % numerical aperture of the objective lens
PRstruct.Lambda = 0.69;                 % center wavelength of the emission band pass filter, unit is micron
PRstruct.RefractiveIndex = 1.52;        % refractive index of the immersion medium
PRstruct.Pupil.phase = zeros(R,R);      % initialize pupil function
PRstruct.Pupil.mag = zeros(R,R);        %
PRstruct.SigmaX = 2;                    % sigmax of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaX), unit is micron
PRstruct.SigmaY = 2;                    % sigmay of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaY), unit is micron
%% generate PSFs
psfobj = PSF_zernike(PRstruct);         % create object from PSF_zernike class
Num = 21;                               % number of PSFs
zpos = linspace(-2,2,Num)';             % z positions of the PSFs , unit is micron
psfobj.Xpos = zeros(Num,1);             % x positions of the PSFs, unit is pixel
psfobj.Ypos = zeros(Num,1);             % y positions of the PSFs, unit is pixel
psfobj.Zpos = zpos;                     % 
psfobj.Boxsize = 32;                    % output size of the PSFs
psfobj.Pixelsize = 0.13;                % pixel size on the sample plane, unit is micron
psfobj.PSFsize = R;                     % image size used for PSF generation
psfobj.nMed = 1.33;                     % refractive index for sample medium

psfobj.precomputeParam();               % generate parameters for Fourier space operation
psfobj.genPupil();                      % generate pupil function 
psfobj.genPSF();                        % generate PSFs
psfobj.scalePSF('normal');              % generate OTF scaled PSFs, 'IMM': PSFs with index mismatch aberration, 'normal': PSFs without index mismatch aberration 
psf = psfobj.ScaledPSFs;                % simulated PSFs
psfobj.genfigs('normal');               % show simulated PSFs, 'IMM': PSFs with index mismatch aberration, 'normal': PSFs without index mismatch aberration
%% calculate CRLB from above PSF model
photon = 8000;                          % total photon count of a single emitter
bg = 20;                                % background photon count of each pixel
crobj = CalCRLB(PRstruct,'zernike');    % create object for CalCRLB class
                                        % input parameters: PRstruct,PSFtype
                                        % PSFtype - 'zernike' : PSFs without index mismatch aberration and generated from PSF_zernike class 
                                        %           'IMM'     : PSFs with index mismatch aberration and generated from PSF_zernike class 
                                        %           'DH'      : double helix PSFs   
crobj.Photon = photon.*ones(Num,1);     %
crobj.Bg = zeros(Num,1)+bg;             %
crobj.Xpos = psfobj.Xpos;               % --------copy parameters from psfobj to crobj-----
crobj.Ypos = psfobj.Ypos;               %
crobj.Zpos = psfobj.Zpos;               %
crobj.Pixelsize = psfobj.Pixelsize;     % 
crobj.PSFobj.PSFsize = psfobj.PSFsize;  %
crobj.PSFobj.nMed = psfobj.nMed;        %
crobj.Boxsize = psfobj.Boxsize;         % -------------------------------------------------
crobj.Deltax = 0.1;                     % increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
crobj.Deltaz = 0.01;                    % increment in z direction for calculating first and second derivative of the objective function, unit is micron 
 
crobj.prepInputparam();                 % generate parameters of PSFs used for CRLB calculation
crobj.calcrlb();                        % calculate CRLB of simulated emitters, given a PSF model
crobj.genfigs();                        % generate plots of theoretical localization precision in x, y and z at z positions defined by 'crobj.Zpos'

% output: crobj.X_STD       - precision in x estimation
%         crobj.Y_STD       - precision in y estimation
%         crobj.Z_STD       - precision in z estimation

