
% (C) Copyright 2018                
%     All rights reserved           
%
% Author: Sheng Liu, May 2018

classdef DipoleField < handle
    % DipoleField class for calculating dipole emission. 
    %   The theoretical calculation is reference to document
    %   'ISM simulation and theory.pdf'
    %
    %   create object: obj = DipoleField();
    %   
    % DipoleField Properties (Input):
    %   NA - 
    %   Lambda - 
    %   nImm - 
    %   nMed - 
    %   nCov - 
    %   Pixelsize - 
    %   PSFsize - 
    %   Zpos - 
    %   x0 - 
    %   y0 - 
    %   DipoleP - 
    %   DipoleAlpha - 
    %   DipoleBeta - 
    %   Exc_x - 
    %   Exc_y - 
    %   Exc_z - 
    %   Nv - 
    %   NV1 - 
    %
    % DipoleField Properties (Output):
    %   Ex - 
    %   Ey - 
    %   I - 
    %
    % DipoleField Methods:
    %   precomputeParam - calculate all parameters defined in the second part of properties
    %   fieldxy - calculate the x and y component of the electric field of dipole emission on the image plane
    %   singledipole - calculate the electric field of a fixed dipole given its dipole moment. 
    %   isodipole - calculate the electric field of two types of rotating dipole.
    %   calDipoleV - calculate x, y, and z components of dipole moment.
    
    properties
        NA;% numerical aperture of objective lens
        Lambda;% emission wavelength
        nImm; % refractive index of immersion medium
        nMed; % refractive index of sample medium
        nCov; % refractive index of cover glass
        Pixelsize; % pixel size at sample plane, in micron
        PSFsize; % original PSF size, in pixels
        Boxsize; % out put PSF size
        Pixelsizefine;% Pixelsize of original PSF
        Zpos; % z positions of cover glass, in micron, could be a scalar or a vector
        Xpos; % x position of the dipole emitter, in pixels, scalar or vector
        Ypos; % y position of the dipole emitter, in pixels, scalar or vector
        Alphas; % polar angle of the dipole emitter, scalar or vector
        Betas; % azimuthal angle of the dipole emitter, scalar or vector
        z0=0; % z position of the dipole emitter, in micron, must be a scalar
        x0=0; % x position of the dipole emitter, in micron, must be a scalar
        y0=0; % y position of the dipole emitter, in micron, must be a scalar
        DipoleP=1; % magnitude of dipole moment
        DipoleAlpha=0; % polar angle of dipole moment
        DipoleBeta=0; % azimuthal angle of dipole moment
        Exc_x; % x component of excitation field
        Exc_y; % y component of excitation field
        Exc_z; % z component of excitation field
        NV=16; % summation number in theta integration
        % NV1 - for summation over istropical dipole distribution, NV1 of
        % sample points of azimuthal angle Beta from 0 to 2*pi, 2*NV1 of
        % sample points of polar angle Alpha from 0 to pi.
        NV1=19; 
    end
    
    properties (SetAccess = private, GetAccess = private)
        % precompute parameters
        Px; % x component of dipole moment
        Py; % y component of dipole moment
        Pz; % z component of dipole moment
        Cos1; % cos(phi_m), a vector of PSFsize x PSFsize elements
        Sin1; % sin(phi_m), a vector of PSFsize x PSFsize elements
        Cos2; % cos(2*phi_m), a vector of PSFsize x PSFsize elements
        Sin2; % sin(2*phi_m), a vector of PSFsize x PSFsize elements
        B0; % J_0(kr), a PSFsize^2 x NV matrix  
        B1; % J_1(kr), a PSFsize^2 x NV matrix  
        B2; % J_2(kr), a PSFsize^2 x NV matrix  
        C0; % a vector of NV elements
        C1; % a vector of NV elements
        C2; % a vector of NV elements
        % Tp - total transmission coefficient of p-polarization (parallel
        % to incident plane) for 3-layer system, from sample medium,
        % coverglass, to immersion medium, a vector of NV elements
        Tp; 
        % Ts - total transmission coefficient of s-polarization
        % (perpendicular to incident plane) for 3-layer system, from sample
        % medium, coverglass, to immersion medium, a vector of NV elements
        Ts; 
        Taup; % a vector of NV elements
        Taus; % a vector of NV elements
        Rous; % a vector of NV elements
        Roup; % a vector of NV elements
        Kz3; % z component of k vector in immersion medium
    end
    
    properties (SetAccess = private, GetAccess = public)
        Ex; % x component of electric field of dipole emission, its a matrix of PSFsize x PSFsize x number of z positions
        Ey; % y component of electric field of dipole emission, its a matrix of PSFsize x PSFsize x number of z positions
        I; % intensity distribution of dipole emission, its a matrix of PSFsize x PSFsize x number of z positions
        PSFs; % point spread function of dipole emitter
    end
    
    methods
        function obj=DipoleField
        end
        
        function precomputeParam(obj)
            % precomputeParam - calculate all parameters defined in the
            % second part of properties.
            %   They are only depend on the optical system, not depend on the
            %   dipole emitter or the excitation field. The naming convention
            %   for multilayer system is : 1 in water, 2 cover glass, 3
            %   immersion medium

            rg=-floor(obj.PSFsize/2)+(1-mod(obj.PSFsize,2))*0.5:floor(obj.PSFsize/2)-(1-mod(obj.PSFsize,2))*0.5;
            [X,Y]=meshgrid(rg,rg);
            Zo=sqrt(X.^2+Y.^2);
            r0=sqrt(obj.x0^2+obj.y0^2);
            phi0=atan2(-obj.y0,-obj.x0);
            k=2*pi*obj.nImm/obj.Lambda;
            phi=atan2(Y,X);
            r=Zo.*obj.Pixelsizefine;% real space
            AngleMax=asin(obj.NA/obj.nImm);
            theta=linspace(0.001,AngleMax,obj.NV);
            
            rshift=sqrt(r.^2+r0^2+2.*r.*r0.*cos(phi-phi0));
            vecr=reshape(rshift,obj.PSFsize^2,1);
            phishift=atan2(r.*sin(phi)+r0*sin(phi0),r.*cos(phi)+r0*cos(phi0));
            vecphi=reshape(phishift,obj.PSFsize^2,1);

            sin_theta3=sin(theta);
            sin_theta1=obj.nImm./obj.nMed.*sin_theta3;
            sin_theta2=obj.nImm./obj.nCov.*sin_theta3;
            
            cos_theta1=sqrt(1-sin_theta1.^2);
            cos_theta2=sqrt(1-sin_theta2.^2);
            cos_theta3=sqrt(1-sin_theta3.^2);
            
            obj.Cos2=cos(2.*vecphi);
            obj.Sin2=sin(2.*vecphi);
            obj.Cos1=cos(vecphi);
            obj.Sin1=sin(vecphi);
            
            obj.B0=[];
            obj.B1=[];
            obj.B2=[];
            for ii=1:obj.NV
                R0=vecr.*k.*sin_theta3(ii);
                obj.B0=cat(2,obj.B0,besselj(0,R0));
                obj.B1=cat(2,obj.B1,besselj(1,R0));
                obj.B2=cat(2,obj.B2,besselj(2,R0));
            end
            % Fresnel coefficient, p parallel, s senkrecht, perpendicular
            
            rou12s=(obj.nMed.*cos_theta1-obj.nCov.*cos_theta2)./(obj.nMed.*cos_theta1+obj.nCov.*cos_theta2);
            tau12s=2*obj.nMed.*cos_theta1./(obj.nMed.*cos_theta1+obj.nCov.*cos_theta2);
            rou23s=(obj.nCov.*cos_theta2-obj.nImm.*cos_theta3)./(obj.nCov.*cos_theta2+obj.nImm.*cos_theta3);
            tau23s=2*obj.nCov.*cos_theta2./(obj.nCov.*cos_theta2+obj.nImm.*cos_theta3);
            
            rou12p=(obj.nCov.*cos_theta1-obj.nMed.*cos_theta2)./(obj.nCov.*cos_theta1+obj.nMed.*cos_theta2);
            tau12p=2*obj.nMed.*cos_theta1./(obj.nCov.*cos_theta1+obj.nMed.*cos_theta2);
            rou23p=(obj.nImm.*cos_theta2-obj.nCov.*cos_theta3)./(obj.nImm.*cos_theta2+obj.nCov.*cos_theta3);
            tau23p=2*obj.nCov.*cos_theta2./(obj.nImm.*cos_theta2+obj.nCov.*cos_theta3);
            
            T0=sqrt(cos_theta3.*obj.nImm./cos_theta1./obj.nMed);
            %T0=1;
            obj.Taup=tau12p.*tau23p.*T0;
            obj.Roup=rou12p.*rou23p;
            obj.Taus=tau12s.*tau23s.*T0;
            obj.Rous=rou12s.*rou23s;
            
            obj.C0=sin_theta3.*cos_theta3./cos_theta1;
            obj.C1=cos_theta1.*obj.C0;
            obj.C2=sin_theta1.*obj.C0;
            obj.Kz3=k.*cos_theta3;
        end
        
        function [Ex,Ey,I]=fieldxy(obj)
            % fieldxy - calculate the x and y component of the electric
            % field of dipole emission on the image plane. 
            %   The subcription '_pxx' means:
            %   p - parallel;
            %   x - x component of dipole moment; 
            %   x - x component of electric field.
%            N = obj.NV;
%             vecE_x = zeros(obj.PSFsize^2,N);
%             vecE_y = zeros(obj.PSFsize^2,N);    
            vecE_x = 0;
            vecE_y = 0;
            for ii=1:obj.NV
                
                vecE_pxx=obj.Px.*obj.C1(ii).*pi.*(obj.B0(:,ii)-obj.Cos2.*obj.B2(:,ii)).*obj.Tp(ii);
                vecE_pyx=-obj.Py.*obj.C1(ii).*pi.*obj.Sin2.*obj.B2(:,ii).*obj.Tp(ii);
                vecE_pzx=-obj.Pz.*obj.C2(ii).*1i.*2.*pi.*obj.Cos1.*obj.B1(:,ii).*obj.Tp(ii);
                vecE_sxx=obj.Px.*obj.C0(ii).*pi.*(obj.B0(:,ii)+obj.Cos2.*obj.B2(:,ii)).*obj.Ts(ii);
                vecE_syx=obj.Py.*obj.C0(ii).*pi.*obj.Sin2.*obj.B2(:,ii).*obj.Ts(ii);
                
                vecE_pxy=-obj.Px.*obj.C1(ii).*pi.*obj.Sin2.*obj.B2(:,ii).*obj.Tp(ii);
                vecE_pyy=obj.Py.*obj.C1(ii).*pi.*(obj.B0(:,ii)+obj.Cos2.*obj.B2(:,ii)).*obj.Tp(ii);
                vecE_pzy=-obj.Pz.*obj.C2(ii).*2.*1i.*pi.*obj.Sin1.*obj.B1(:,ii).*obj.Tp(ii);
                vecE_sxy=obj.Px.*obj.C0(ii).*pi.*obj.Sin2.*obj.B2(:,ii).*obj.Ts(ii);
                vecE_syy=obj.Py.*obj.C0(ii).*pi.*(obj.B0(:,ii)-obj.Cos2.*obj.B2(:,ii)).*obj.Ts(ii);
                                
%                 vecE_x(:,ii) = vecE_pxx+vecE_pyx+vecE_pzx+vecE_sxx+vecE_syx;
%                 vecE_y(:,ii) = vecE_pxy+vecE_pyy+vecE_pzy+vecE_sxy+vecE_syy;
                vecE_x = vecE_x + vecE_pxx+vecE_pyx+vecE_pzx+vecE_sxx+vecE_syx;
                vecE_y = vecE_y + vecE_pxy+vecE_pyy+vecE_pzy+vecE_sxy+vecE_syy;
            end
%             vecE_x = sum(vecE_x,2);
%             vecE_y = sum(vecE_y,2);
            Ex=reshape(vecE_x,obj.PSFsize,obj.PSFsize);
            Ey=reshape(vecE_y,obj.PSFsize,obj.PSFsize);
            I=Ex.*conj(Ex)+Ey.*conj(Ey);
            
        end
        function singledipole(obj)
            % singledipole - calculate the electric field of a fixed dipole
            % given its dipole moment.
            %   If 'obj.Zpos' is a vector, it will loop over all z positions
            %   and calculate the electric field at each z position. The
            %   output results 'Ex', 'Ey' and 'I' are 3D matrices
            %   of PSFsize x PSFsize x number of z positions
            %
            %   see also fieldxy
            Ex0=[];
            Ey0=[];
            I0=[];
            for ii=1:length(obj.z0)
                z=obj.z0(ii);% um
                beta3=obj.Kz3.*z;
                obj.Tp=obj.Taup./(obj.Roup.*exp(-1i.*beta3)+exp(-1i.*beta3));
                obj.Ts=obj.Taus./(obj.Rous.*exp(-1i.*beta3)+exp(-1i.*beta3));
                [Exi,Eyi,Ii]=obj.fieldxy();
                Ex0=cat(3,Ex0,Exi);
                Ey0=cat(3,Ey0,Eyi);
                I0=cat(3,I0,Ii);
            end
            obj.Ex=Ex0;
            obj.Ey=Ey0;
            obj.I=I0;
        end
        
        function genPSF(obj,bin)
            N = length(obj.Xpos);
            obj.PSFs = zeros(obj.Boxsize,obj.Boxsize,N);
            obj.Pixelsizefine = obj.Pixelsize/bin;
            obj.PSFsize = obj.Boxsize*bin;
           for ii = 1:N
                obj.z0 = obj.Zpos(ii);% um
                obj.x0 = obj.Xpos(ii)*obj.Pixelsize;% um
                obj.y0 = obj.Ypos(ii)*obj.Pixelsize;% um
                obj.DipoleAlpha = obj.Alphas(ii); % polar angle
                obj.DipoleBeta = obj.Betas(ii);% azimuthal angle
                
                obj.precomputeParam();
                obj.calDipoleV();
                obj.singledipole();
                psf = binimage(obj.I,bin);
                obj.PSFs(:,:,ii) = psf./sum(psf(:));
            end
        end
        
        function genPSF1(obj,bin,psfobj)
            N = length(obj.Xpos);
            psfs = zeros(obj.Boxsize,obj.Boxsize,N);
            parfor ii = 1:N
                tmpobj = psfobj;
                tmpobj.Pixelsizefine = obj.Pixelsize/bin;
                tmpobj.PSFsize = obj.Boxsize*bin;
                
                tmpobj.z0 = obj.Zpos(ii);% um
                tmpobj.x0 = obj.Xpos(ii)*obj.Pixelsize;% um
                tmpobj.y0 = obj.Ypos(ii)*obj.Pixelsize;% um
                tmpobj.DipoleAlpha = obj.Alphas(ii); % polar angle
                tmpobj.DipoleBeta = obj.Betas(ii);% azimuthal angle
                
                tmpobj.precomputeParam();
                tmpobj.calDipoleV();
                tmpobj.singledipole();
                psf = binimage(tmpobj.I,bin);
                psfs(:,:,ii) = psf./sum(psf(:));
            end
            obj.PSFs = psfs;
        end
        function isodipole(obj,dipoleType)
            % isodipole - calculate the electric field of two types of
            % rotating dipole.
            %   One is 'fast' rotating dipole, which has no memory of the
            %   polarization direction of the electric field, therefore the
            %   dipole moment only scaled by the magnitude of the electric
            %   field, for simplicity, the scale is set to 1. 
            %
            %   One is 'slow' rotating dipole, which dipole moment is depend
            %   both on the polarization and the magnitude of the electric
            %   field. 
            %
            %   The resulting intensity distribution is an sum over all the
            %   sample points of dipole orientation in 4pi degree. Because
            %   of the sum, the x and y component of electric field is
            %   meaningless, therefore set 'Ex' and 'Ey' empty.
            %
            %   Input parameter: dipoleType - 'fast' or 'slow'
            %
            %   see also singledipole
            a=linspace(0,pi,obj.NV1);
            b=linspace(0,2*pi,2*obj.NV1);
            [aa,bb]=meshgrid(a,b);
            aa=aa(:);
            bb=bb(:);
            I0=zeros(obj.PSFsize,obj.PSFsize,length(obj.z0));
            switch dipoleType
                case 'fast'
                    c=[0;1];
                case 'slow'
                    c=[1;0];
            end
            for ii=1:numel(aa)
                px0=sin(aa(ii))*cos(bb(ii));
                py0=sin(aa(ii))*sin(bb(ii));
                pz0=cos(aa(ii));
                G=[px0,py0,pz0]*[obj.Exc_x,obj.Exc_y,obj.Exc_z]';
                G=[G,1]*c;
                obj.Px=px0*G;
                obj.Py=py0*G;
                obj.Pz=pz0*G;
                obj.singledipole();
                I0=I0+obj.I;
            end
            obj.I=I0;
            obj.Ex=[];
            obj.Ey=[];
        end
        
        function calDipoleV(obj)
            % calDipoleV - calculate x, y, and z components of dipole
            % moment.
            %   This is only used in calculating the electric field of a
            %   single dipole.
            %
            %   see also singledipole
            P=obj.DipoleP; 
            a=obj.DipoleAlpha*pi/180; 
            b=obj.DipoleBeta*pi/180; 
            obj.Px=P*sin(a)*cos(b);
            obj.Py=P*sin(a)*sin(b);
            obj.Pz=P*cos(a);

        end
    end
    
end

