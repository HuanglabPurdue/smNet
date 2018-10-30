
% (C) Copyright 2018                
%     All rights reserved           
%
% Author: Sheng Liu, May 2018

classdef CalCRLB < handle
    % CalCRLB class for calculating CRLB of simulated emitters, given a PSF model defined by PRstruct
    %   create object: obj = CalCRLB(PRstruct,PSFtype)
    %   inputs required: PRstruct - define necessary parameters for a PSF model
    %                   PSFtype - type of method to generate PSFs for CRLB
    %                   calculation, options are 'pupil' and 'zernike'
    %
    % CalCRLB Properties (Implemented object):
    %   PSFobj - 
    %
    % CalCRLB Properties (Input):
    %   Bg - 
    %   Boxsize - 
    %   Deltax - 
    %   Deltaz - 
    %   PRstruct - 
    %   Photon - 
    %   Pixelsize - 
    %   Xpos - 
    %   Ypos - 
    %   Zpos - 
    %
    % CalCRLB Properties (Output):
    %   CRLB - 
    %   X_STD - 
    %   Y_STD - 
    %   Z_STD - 
    %   Photon_STD - 
    %   Bg_STD - 
    %
    % CalCRLB Methods:
    %   prepInputparam - generate parameters of PSFs used for CRLB calculation
    %   calcrlb - calculate CRLB of simulated emitters, given a PSF model
    %   genfigs - generate plots of theoretical localization precision in x, y and z at z positions defined by 'Zpos'
    %
    %   see also PSF_pupil
    properties
        % PRstruct - define necessary parameters for a PSF model
        %   NA
        %   Lambda
        %   RefractiveIndex
        %   Pupil: phase retrieved pupil function
        %           phase: phase image
        %           mag: magnitude image
        %   Zernike_phase: coefficient of zernike polynomials representing the pupil phase
        %   Zernike_mag: coefficient of zernike polynomials representing the pupil phase
        %   SigmaX: sigmax of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaX), unit is micron
        %   SigmaY: sigmay of Gaussian filter for OTF rescale, unit is 1/micron in k space, the conversion to real space is 1/(2*pi*SigmaY), unit is micron
        PRstruct;
        PSFobj;% object of PSF_pupil or PSF_zernike class, used for generating PSFs 
        Xpos;% x positions of simulated emitters, a vector of N elements, unit is pixel
        Ypos;% y positions of simulated emitters, a vector of N elements, unit is pixel
        Zpos;% z positions of simulated emitters, a vector of N elements, unit is micron
        Photon;% photon counts of simulated emitters, a vector of N elements
        Bg;% % background photon counts of simulated emitters, a vector of N elements
        Pixelsize;% pixel size at sample plane, unit is micron
        Boxsize;% image size of simulated emitter
        PSFtype;% type of method to generate PSFs for CRLB
        Deltax;% increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
        Deltaz;% increment in z directions forcalculating first and second derivative of the objective function, unit is micron 
        PN=5; % number of fitting parameters in CRLB calculation
    end
    properties (SetAccess = private, GetAccess = public)
        Xin; % parameters of PSFs, a N*PN x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
        PSFI;
    end
    % output parameters
    properties (SetAccess = private, GetAccess = public)
        CRLB; % CRLB of simulated emmiters, a N x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
        X_STD; % theoretical localization precision in X dimension, a vector of N elements, unit is pixel
        Y_STD; % theoretical localization precision in Y dimension, a vector of N elements, unit is pixel
        Z_STD; % theoretical localization precision in Z dimension, a vector of N elements, unit is micron
        Photon_STD; % theoretical localization precision in photon count, a vector of N elements
        Bg_STD; % theoretical localization precision in background count, a vector of N elements
    end
    
    methods
        function obj = CalCRLB(PRstruct,PSFtype)
            switch PSFtype
                case 'zernike'
                    obj.PSFobj = PSF_zernike(PRstruct);
                case 'IMM'
                    obj.PSFobj = PSF_zernike(PRstruct);
                case 'DH'
                    obj.PSFobj = PSF_DH(PRstruct);
            end
            obj.PSFtype = PSFtype;
        end
        
        function prepInputparam(obj)
            % prepInputparam - generate parameters of PSFs use for CRLB calculation.
            %   the out put is Xin. it is a N*PN x PN matrix, N is the
            %   number of elements in Xpos. PN is the number of fitting
            %   parameters, including x, y, z, photon, bg
            N=numel(obj.Xpos);
            pN=obj.PN;
            PIx0=zeros(1,pN);
            PIy0=zeros(1,pN);
            PIz0=zeros(1,pN);
            PIi0=zeros(1,pN);
            PIi=ones(1,pN);
            
            PIx0(3)=obj.Deltax;
            PIy0(4)=obj.Deltax;
            PIz0(5)=obj.Deltaz;
            PIi0(2)=1;
            PIi(2)=0;
            
            x=[];
            y=[];
            I=[];
            bg=[];
            z=[];
            x0=cat(2,obj.Xpos,obj.Ypos,obj.Photon,obj.Bg,obj.Zpos);
            for t=1:pN
                x = cat(2,x,x0(:,1)+PIx0(t));
                y = cat(2,y,x0(:,2)+PIy0(t));
                I = cat(2,I,x0(:,3).*PIi(t)+PIi0(t));
                bg = cat(2,bg,x0(:,4).*PIi(t));
                z = cat(2,z,x0(:,5)+PIz0(t));
            end
            obj.Xin=cat(2,reshape(x',N*pN,1),reshape(y',N*pN,1),reshape(I',N*pN,1),reshape(bg',N*pN,1),reshape(z',N*pN,1));
            
        end
        
        function calcrlb(obj)
            % calcrlb - calculate CRLB of simulated emitters, given a PSF model.
            %   It uses PSF_pupil or PSF_zernike class to generate PSFs
            %   from given parameters in 'Xin'
            %
            %   see also PSF_pupil
            obj.PSFobj.Xpos=obj.Xin(:,1);
            obj.PSFobj.Ypos=obj.Xin(:,2);
            if strcmp(obj.PSFtype,'IMM')
                obj.PSFobj.ZposMed = obj.Xin(:,5);
            else
                obj.PSFobj.Zpos=obj.Xin(:,5);
            end
            obj.PSFobj.Boxsize=obj.Boxsize;
            obj.PSFobj.Pixelsize=obj.Pixelsize;
            if ~sum(strcmp(obj.PSFtype,{'interp','gauss'}))
                obj.PSFobj.precomputeParam();
            end
            switch obj.PSFtype
                case 'zernike'
                    obj.PSFobj.genPupil();
                    obj.PSFobj.genPSF();
                    obj.PSFobj.scalePSF('normal');
                case 'IMM'
                    obj.PSFobj.genPupil();
                    obj.PSFobj.genIMMPSF();
                    obj.PSFobj.scalePSF('IMM');
                case 'DH'
                    obj.PSFobj.genPupil();
                    obj.PSFobj.genPSF();
            end
            % calculate Fisher Information matrix
            N = numel(obj.Xpos);
            pN = obj.PN;
            funFi = zeros(obj.Boxsize,obj.Boxsize,pN);
            FisherM = zeros(pN,pN);
            xVar = zeros(N,obj.PN);
            psfI = obj.PSFobj.ScaledPSFs;
            obj.PSFI = psfI;
            tmp = obj.PSFobj.PRstruct.Pupil.mag;
            normf = sum(sum(tmp.^2,1),2);
            switch obj.PSFtype
                case 'zernike'
                    psf = psfI./normf;
                case 'IMM'
                    psf = psfI./normf;
                case 'DH'
                    psf = psfI;
            end
            
            for s=0:N-1
                t=s+1;
                %x
                funFi(:,:,1)=obj.Photon(t).*(psf(:,:,s*pN+3)-psf(:,:,s*pN+1))./obj.Deltax;
                %y
                funFi(:,:,2)=obj.Photon(t).*(psf(:,:,s*pN+4)-psf(:,:,s*pN+1))./obj.Deltax;
                %I
                funFi(:,:,3)=psf(:,:,s*pN+2);
                %bg
                funFi(:,:,4)=1;
                %z
                funFi(:,:,5)=obj.Photon(t).*(psf(:,:,s*pN+5)-psf(:,:,s*pN+1))./obj.Deltaz;
                
                for j=1:pN
                    for k=1:pN
                        psfIni=psf(:,:,s*pN+1).*obj.Photon(t)+obj.Bg(t);
                        FisherM(j,k)=sum(sum(funFi(:,:,j).*funFi(:,:,k)./psfIni));
                    end
                end
                LowerBi=inv(FisherM);
                xVar(t,:)=diag(LowerBi)';
            end
            obj.CRLB=abs(xVar);
            obj.X_STD=sqrt(obj.CRLB(:,1));
            obj.Y_STD=sqrt(obj.CRLB(:,2));
            obj.Z_STD=sqrt(obj.CRLB(:,5));
            obj.Photon_STD=sqrt(obj.CRLB(:,3));
            obj.Bg_STD=sqrt(obj.CRLB(:,4));
        end
        
        function genfigs(obj)
            % genfigs - generate plots of theoretical localization
            % precision in x, y and z at z positions defined by Zpos
            figure('position',[100,200,500,300])
            plot(obj.Zpos.*1e3,obj.X_STD.*obj.Pixelsize.*1e3,'r.-')
            hold on
            plot(obj.Zpos.*1e3,obj.Y_STD.*obj.Pixelsize.*1e3,'b.-')
            plot(obj.Zpos.*1e3,obj.Z_STD.*1e3,'g.-')
            axis tight;
            xlabel('z positions (nm)')
            ylabel('precison from CRLB (nm)')
            legend('\sigmax','\sigmay','\sigmaz')
        end
    end
    
end
