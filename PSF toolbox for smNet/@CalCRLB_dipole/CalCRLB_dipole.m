
% (C) Copyright 2018                
%     All rights reserved           
%
% Author: Sheng Liu, May 2018

classdef CalCRLB_dipole < handle
    
    properties
        PSFobj;% object of PSF_pupil or PSF_zernike class, used for generating PSFs
        Xpos;% x positions of simulated emitters, a vector of N elements, unit is pixel
        Ypos;% y positions of simulated emitters, a vector of N elements, unit is pixel
        Zpos;% z positions of simulated emitters, a vector of N elements, unit is micron
        Alphas; % polar angle of the dipole emitter
        Betas;% azimuthal angle of the dipole emitter
        Photon;% photon counts of simulated emitters, a vector of N elements
        Bg;% background photon counts of simulated emitters, a vector of N elements
        Pixelsize;% pixel size at sample plane, unit is micron
        Boxsize;% image size of simulated emitter
        PSFtype;% type of method to generate PSFs for CRLB
        Deltax;% increment in x and y directions for calculating first and second derivative of the objective function, unit is pixel
        Deltaz;% increment in z directions forcalculating first and second derivative of the objective function, unit is micron
        Deltaa;% increment in polar angle
        Deltab;% increment in azimuthal angle
        PN = 6;
        Bin = 2;
    end
    
    properties (SetAccess = private, GetAccess = public)
        Xin; % parameters of PSFs, a N*PN x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
    end
    % output parameters
    properties (SetAccess = private, GetAccess = public)
        CRLB; % CRLB of simulated emmiters, a N x PN matrix, N is the number of elements in Xpos. PN is the number of fitting parameters, including x, y, z, photon, bg  
        X_STD; % theoretical localization precision in X dimension, a vector of N elements, unit is pixel
        Y_STD; % theoretical localization precision in Y dimension, a vector of N elements, unit is pixel
        Z_STD; % theoretical localization precision in Z dimension, a vector of N elements, unit is micron
        Photon_STD; % theoretical localization precision in photon count, a vector of N elements
        Bg_STD; % theoretical localization precision in background count, a vector of N elements
        Alpha_STD; % theoretical localization precision in polar angle, a vector of N elements
        Beta_STD; % theoretical localization precision in azimuthal angle, a vector of N elements
    end
    
    methods
        function obj = CalCRLB_dipole()
            obj.PSFobj = DipoleField();
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
            PIa0=zeros(1,pN);
            PIb0=zeros(1,pN);
            
            PIx0(2)=obj.Deltax;
            PIy0(3)=obj.Deltax;
            PIz0(4)=obj.Deltaz;
            PIa0(5)=obj.Deltaa;
            PIb0(6)=obj.Deltab;
            
            x=[];
            y=[];
            z=[];
            a=[];
            b=[];
            x0=cat(2,obj.Xpos,obj.Ypos,obj.Zpos,obj.Alphas,obj.Betas);
            for t=1:pN
                x = cat(2,x,x0(:,1)+PIx0(t));
                y = cat(2,y,x0(:,2)+PIy0(t));
                z = cat(2,z,x0(:,3)+PIz0(t));
                a = cat(2,a,x0(:,4)+PIa0(t));
                b = cat(2,b,x0(:,5)+PIb0(t));
            end
            obj.Xin=cat(2,reshape(x',N*pN,1),reshape(y',N*pN,1),reshape(z',N*pN,1),reshape(a',N*pN,1),reshape(b',N*pN,1));

        end
        
        function calcrlb(obj)
            % calcrlb - calculate CRLB of simulated emitters, given a PSF model.
            %   It uses PSF_pupil or PSF_zernike class to generate PSFs
            %   from given parameters in 'Xin'
            %
            %   see also PSF_pupil
            obj.PSFobj.Xpos = obj.Xin(:,1);
            obj.PSFobj.Ypos = obj.Xin(:,2);
            obj.PSFobj.Zpos = obj.Xin(:,3);
            obj.PSFobj.Alphas = obj.Xin(:,4);
            obj.PSFobj.Betas = obj.Xin(:,5);
            obj.PSFobj.Boxsize = obj.Boxsize;
            obj.PSFobj.Pixelsize = obj.Pixelsize;
            obj.PSFobj.genPSF(obj.Bin);
            
            % calculate Fisher Information matrix
            N=numel(obj.Xpos);
            pN=obj.PN;
            pN0 = pN+1;
            funFi=zeros(obj.Boxsize,obj.Boxsize,pN0);
            FisherM=zeros(pN0,pN0);
            xVar=zeros(N,pN0);
            psf=obj.PSFobj.PSFs;
            for s=0:N-1
                t=s+1;
                    %x
                    funFi(:,:,1)=obj.Photon(t).*(psf(:,:,s*pN+2)-psf(:,:,s*pN+1))./obj.Deltax;
                    %y
                    funFi(:,:,2)=obj.Photon(t).*(psf(:,:,s*pN+3)-psf(:,:,s*pN+1))./obj.Deltax;
                    %z
                    funFi(:,:,3)=obj.Photon(t).*(psf(:,:,s*pN+4)-psf(:,:,s*pN+1))./obj.Deltaz;
                    %a
                    funFi(:,:,4)=obj.Photon(t).*(psf(:,:,s*pN+5)-psf(:,:,s*pN+1))./obj.Deltaa;
                    %z
                    funFi(:,:,5)=obj.Photon(t).*(psf(:,:,s*pN+6)-psf(:,:,s*pN+1))./obj.Deltab;
                    %I
                    funFi(:,:,6)=psf(:,:,s*pN+1);
                    %bg
                    funFi(:,:,7)=1;
                    
                    for j=1:pN0
                        for k=1:pN0
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
            obj.Z_STD=sqrt(obj.CRLB(:,3));
            obj.Alpha_STD = sqrt(obj.CRLB(:,4));
            obj.Beta_STD = sqrt(obj.CRLB(:,5));
            obj.Photon_STD=sqrt(obj.CRLB(:,6));
            obj.Bg_STD=sqrt(obj.CRLB(:,7));
        end
        
        function genfigs(obj)
            % genfigs - generate plots of theoretical localization
            % precision in x, y and z at z positions defined by Zpos
            figure('position',[100,200,500*3,300*2])
            subplot(231)
            plot(obj.Zpos.*1e3,obj.X_STD.*obj.Pixelsize.*1e3,'r.-')
            hold on
            plot(obj.Zpos.*1e3,obj.Y_STD.*obj.Pixelsize.*1e3,'b.-')
            plot(obj.Zpos.*1e3,obj.Z_STD.*1e3,'g.-')
            axis tight;
            xlabel('z positions (nm)')
            ylabel('precison from CRLB (nm)')
            legend('\sigmax','\sigmay','\sigmaz')
            
            subplot(232)
            plot(obj.Alphas,obj.X_STD.*obj.Pixelsize.*1e3,'r.-')
            hold on
            plot(obj.Alphas,obj.Y_STD.*obj.Pixelsize.*1e3,'b.-')
            plot(obj.Alphas,obj.Z_STD.*1e3,'g.-')
            axis tight;
            xlabel('polar angle (radians)')
            ylabel('precison from CRLB (nm)')
            legend('\sigmax','\sigmay','\sigmaz')
            
            subplot(233)
            plot(obj.Betas,obj.X_STD.*obj.Pixelsize.*1e3,'r.-')
            hold on
            plot(obj.Betas,obj.Y_STD.*obj.Pixelsize.*1e3,'b.-')
            plot(obj.Betas,obj.Z_STD.*1e3,'g.-')
            axis tight;
            xlabel('azimutha angle (radians)')
            ylabel('precison from CRLB (nm)')
            legend('\sigmax','\sigmay','\sigmaz')
            
            subplot(234)
            plot(obj.Zpos.*1e3,obj.Alpha_STD,'r.-')
            hold on
            plot(obj.Zpos.*1e3,obj.Beta_STD,'b.-')
            axis tight;
            xlabel('z positions (nm)')
            ylabel('precison from CRLB (nm)')
            legend('\sigma\alpha','\sigma\beta')

            subplot(235)
            plot(obj.Alphas,obj.Alpha_STD,'r.-')
            hold on
            plot(obj.Alphas,obj.Beta_STD,'b.-')
            axis tight;
            xlabel('polar angle (radians)')
            ylabel('precison from CRLB (nm)')
            legend('\sigma\alpha','\sigma\beta')
            
            subplot(236)
            plot(obj.Betas,obj.Alpha_STD,'r.-')
            hold on
            plot(obj.Betas,obj.Beta_STD,'b.-')
            axis tight;
            xlabel('azimuthal angle (radians)')
            ylabel('precison from CRLB (nm)')
            legend('\sigma\alpha','\sigma\beta')
        end
    end
    
end






