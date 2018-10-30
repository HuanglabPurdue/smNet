
% (C) Copyright 2018                
%     All rights reserved           
%
% Author: Sheng Liu, May 2018

classdef PSF_DH < handle
    % PSF_DH class for generating double-helix PSFs
    %   create object: obj = PSF_DH(PRstruct)
    %   input required: PRstruct - define necessary parameters for a PSF model
    
    properties
        % PRstruct - define necessary parameters for a PSF model
        %   NA
        %   Lambda
        %   RefractiveIndex
        %   Pupil: pupil function
        %           phase: phase image
        %           mag: magnitude image
        PRstruct;
        Xpos;% x positions of simulated emitters, a vector of N elements, unit is pixel
        Ypos;% y positions of simulated emitters, a vector of N elements, unit is pixel
        Zpos;% z positions of simulated emitters relative to the immersion medium, a vector of N elements, unit is micron
        PSFsize; % image size of the SamplePSF
        Boxsize; % image size of out put PSF
        Pixelsize;% pixel size at sample plane, unit is micron
        Magnification; % magnification of the system
        Zlim; % the upper and lower limit of z positions, [zmin,zmax], unit in micron
    end
    
    properties(SetAccess = private, GetAccess = private)
        % precompute images for k space operation
        Zo;% r coordinates of out put PSF, it's a image of PSFsize x PSFsize
        k_r;% k_r coordinates of out put OTF, it's a image of PSFsize x PSFsize
        k_z;% k_z coordinates of out put OTF, it's a image of PSFsize x PSFsize
        Phi;% phi coordinates out put PSF, it's a image of PSFsize x PSFsize
        NA_constrain;% a circle function defining the limit of k_r, it's a image of PSFsize x PSFsize
    end
    
     properties (SetAccess = private, GetAccess = public)
        % ScaledPSFs - out put PSFs from Fourier transform of the pupil function,
        % it's a 3D matrix of Boxsize x Boxsize x N, N is the number of
        % elements in Xpos.
        ScaledPSFs;
        % Pupil - pupil function generated from a set of zernike polynomials
        %           phase: phase image of PSFsize x PSFsize
        %           mag: magnitude image of PSFsize x PSFsize
        Pupil;
     end
     
      methods
        function obj=PSF_DH(PRstruct)
            obj.PRstruct = PRstruct;
        end
        
         function precomputeParam(obj)
            % precomputeParam - generate images for k space operation, and saved in
            % precomputed parameters.
            [X,Y]=meshgrid(-obj.PSFsize/2:obj.PSFsize/2-1,-obj.PSFsize/2:obj.PSFsize/2-1);
            obj.Zo=sqrt(X.^2+Y.^2);
            scale=obj.PSFsize*obj.Pixelsize;
            obj.k_r=obj.Zo./scale;
            obj.Phi=atan2(Y,X);
            n=obj.PRstruct.RefractiveIndex;
            Freq_max=obj.PRstruct.NA/obj.PRstruct.Lambda;
            obj.NA_constrain=obj.k_r<Freq_max;
            obj.k_z=sqrt((n/obj.PRstruct.Lambda)^2-obj.k_r.^2).*obj.NA_constrain;
         end
         
         function genPupil(obj)
             R = obj.PSFsize;
             n0 = obj.PRstruct.RefractiveIndex;
             lambda = obj.PRstruct.Lambda;
             NA = obj.PRstruct.NA;
             pixelsize = obj.Pixelsize*obj.Magnification;
             Magnify = obj.Magnification;
             zRange = linspace(obj.Zlim(1),obj.Zlim(2),21);
             plotflag = 0;
             [pupil_phase,pupil_mag] = genphaseplate(R,n0,lambda,NA,pixelsize,Magnify,zRange,plotflag);
             obj.Pupil.phase=pupil_phase;
             obj.Pupil.mag=pupil_mag;
             obj.PRstruct.Pupil.phase=pupil_phase;
             obj.PRstruct.Pupil.mag=pupil_mag;
             
         end
         
         function genPSF(obj)
            % genPSF - generate PSFs from the given pupil function.
            %   The PSFs are directly calculated from the Fourier transform
            %   of pupil functions modified by shift phase in x, y and
            %   defocus phase in z. The out put is 'PSFs'
            N=numel(obj.Xpos);
            R=obj.PSFsize;
            Ri=obj.Boxsize;
            psfs=zeros(Ri,Ri,N);
            for ii=1:N
                shiftphase=-obj.k_r.*cos(obj.Phi).*obj.Xpos(ii).*obj.Pixelsize-obj.k_r.*sin(obj.Phi).*obj.Ypos(ii).*obj.Pixelsize;
                shiftphaseE=exp(-1i.*2.*pi.*shiftphase);
                defocusphaseE=exp(2.*pi.*1i.*obj.Zpos(ii).*obj.k_z);
                pupil_complex=obj.Pupil.mag.*obj.Pupil.phase.*shiftphaseE.*defocusphaseE;
                psfA=abs(fftshift(fft2(pupil_complex)));
                Fig2=psfA.^2;
                realsize0=floor(Ri/2);
                realsize1=ceil(Ri/2);
                startx=-realsize0+R/2+1;endx=realsize1+R/2;
                starty=-realsize0+R/2+1;endy=realsize1+R/2;
                psfs(:,:,ii)=Fig2(startx:endx,starty:endy)./sum(Fig2(:));
            end
            obj.ScaledPSFs=psfs;
        end
         
        
      end
   
end