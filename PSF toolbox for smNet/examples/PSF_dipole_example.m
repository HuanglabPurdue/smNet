%------ Demo code for simulation of PSFs of fixed dipole emitter------------
% software requirement: Matlab R2015a or later
%                       Dipimage toolbox 2.7 or later
% system requirement:   CPU Intel Core i7
%                       32 GB RAM 
% Note:this example code demonstrates the usage of DipoleField and CalCRLB_dipole class, 
%      type 'help DipoleField' and 'help CalCRLB_dipole' in the command window for more information
%
% (C) Copyright 2018               
%     All rights reserved          
%
% Author: Sheng Liu, May 2018

%% generate PSFs
                                % create parameter space for dipole emitter
[zpos,~] = meshgrid(linspace(-0.5,0.5,5),linspace(10,90,9));
[beta,alpha] = meshgrid(linspace(0,180,5),linspace(10,90,9));
Num = numel(alpha);             % number of PSFs         
psfobj = DipoleField();         % create object from DipoleField class
psfobj.Xpos = zeros(Num,1);     % x positions of the PSFs, unit is pixel
psfobj.Ypos = zeros(Num,1);     % y positions of the PSFs, unit is pixel
psfobj.Zpos = zpos(:);          % z positions of the PSFs, unit is micron
psfobj.Alphas = alpha(:);       % polar angles of the dipole emitters, unit is degree, range is [0,90]
psfobj.Betas = beta(:);         % azimuthal angles of the dipole emitters, unit is degree, range is [0,360]
psfobj.Boxsize = 32;            % output size of the PSFs
psfobj.NV = 64;                 % division number in 1D numerical integration
psfobj.Pixelsize = 0.05;        % pixel size on the sample plane, unit is micron
psfobj.NA = 1.4;                % numerical aperture of the objective lens
psfobj.nImm = 1.52;             % refractive index of the immersion medium
psfobj.nMed = 1.35;             % refractive index of the sample medium
psfobj.nCov = 1.52;             % refractive index of the coverglass
psfobj.Lambda = 0.69;           % center wavelength of the emission band pass filter, unit is micron
psfobj.DipoleP = 1;             % dipole moment, for fixed dipole emitter and circular polarized illumination, this value doesn't matter

bin = 2;                        % upsampling number used for PSF generation, the pixel size will become psfobj.Pixelsize/bin
psfobj.genPSF(bin);             % generate PSFs
psf = psfobj.PSFs;              % simulated PSFs
%% show PSF images 
xN = size(alpha,1);
yN = size(alpha,2);
h = figure;
h.Position = [200,300,100*xN,105*yN];
for n = 1:Num
    if (n<=xN*yN)
        ha = axes;
        ii = mod(n,xN);
        jj = floor(n/xN)+1;
        if ii == 0
            ii = xN;
            jj = n/xN;
        end
        psfi = psfobj.PSFs(:,:,n);
        ha.Position = [(ii-1)/xN,(yN-jj)/yN,1/xN,1/yN];
        image(psfi,'cdatamapping','scale');axis off;axis equal;
        text(2,3, ['\alpha=',num2str(psfobj.Alphas(n),3),'^o',', \beta=',num2str(psfobj.Betas(n),3),'^o'],'color',[1,1,1],'fontsize',12);
        text(2,psfobj.Boxsize-3, ['z=',num2str(psfobj.Zpos(n),3),'\mum'],'color',[1,1,1],'fontsize',12);
    end
end
colormap(jet);
h.InvertHardcopy = 'off';
h.Position = [200,300,105*xN,105*yN];

%% calculate CRLB from above PSF model, this might take several minutes
photon = 4000;                          % total photon count of a single emitter
bg = 5;                                 % background photon count of each pixel
crobj = CalCRLB_dipole();               % create object from CalCRLB_dipole class
crobj.Bin = 2;                          % upsampling number used for PSF generation
crobj.Photon = photon.*ones(Num,1);     %
crobj.Bg = zeros(Num,1)+bg;             %
crobj.Pixelsize = psfobj.Pixelsize;     % --------copy parameters from psfobj to crobj-----
crobj.Xpos = psfobj.Xpos;               %
crobj.Ypos = psfobj.Ypos;               %
crobj.Zpos = psfobj.Zpos;               %
crobj.Alphas = psfobj.Alphas;           %
crobj.Betas = psfobj.Betas;             %
crobj.Boxsize = psfobj.Boxsize;         %
crobj.PSFobj.NV = psfobj.NV;            %
crobj.PSFobj.NA = psfobj.NA;            %
crobj.PSFobj.nImm = psfobj.nImm;        %
crobj.PSFobj.nMed = psfobj.nMed;        %
crobj.PSFobj.nCov = psfobj.nCov;        %
crobj.PSFobj.Lambda = psfobj.Lambda;    % 
crobj.PSFobj.DipoleP = psfobj.DipoleP;  % --------------------------------------------------
crobj.Deltax = 0.1;                     % increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
crobj.Deltaz = 0.01;                    % increment in z direction for calculating first and second derivative of the objective function, unit is micron 
crobj.Deltaa = 1;                       % increment in polar angle direction for calculating first and second derivative of the objective function, unit is degree 
crobj.Deltab = 1;                       % increment in azimuthal angle direction for calculating first and second derivative of the objective function, unit is degree 
 
crobj.prepInputparam();                 % generate parameters of PSFs used for CRLB calculation
crobj.calcrlb();                        % calculate CRLB of simulated emitters, given a PSF model

% output: crobj.X_STD       - precision in x estimation
%         crobj.Y_STD       - precision in y estimation
%         crobj.Z_STD       - precision in z estimation
%         crobj.Alpha_STD   - precision in polar angle estimation
%         crobj.Beta_STD    - precision in azimuthal angle estimation
%% show results of CRLB calculation
h = figure;
h.Position = [100,50,400*3,300*2];

subplot(231)
surf(alpha,zpos.*1e3,reshape(crobj.X_STD.*crobj.Pixelsize.*1e3,xN,yN))
xlabel('polar angle (degree)')
ylabel('z position (nm)')
zlabel('\sigma_x (nm)')
subplot(232)
surf(alpha,zpos.*1e3,reshape(crobj.Y_STD.*crobj.Pixelsize.*1e3,xN,yN))
xlabel('polar angle (degree)')
ylabel('z position (nm)')
zlabel('\sigma_y (nm)')
subplot(233)
surf(alpha,zpos.*1e3,reshape(crobj.Z_STD.*1e3,xN,yN))
xlabel('polar angle (degree)')
ylabel('z position (nm)')
zlabel('\sigma_z (nm)')
subplot(234)
surf(alpha,zpos.*1e3,reshape(crobj.Alpha_STD,xN,yN))
xlabel('polar angle (degree)')
ylabel('z position (nm)')
zlabel('\sigma_\alpha (\^o)')
subplot(235)
surf(alpha,zpos.*1e3,reshape(crobj.Beta_STD,xN,yN))
xlabel('polar angle (degree)')
ylabel('z position (nm)')
zlabel('\sigma_\beta (\^o)')

