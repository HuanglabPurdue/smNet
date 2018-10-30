% generate gauss-laguerre mode in z stack
function [U]=GLmode(m,n,R,n0,lambda,NA,pixelsize,Magnify,z)
Diffsize=1.22*lambda/NA;% diffraction limited spot size
w0=0.5*Diffsize; %beam waist of GL00 mode
% theta0=asin(NA/n0);
% w0=lambda/n0/pi/tan(theta0);
z0=pi*n0*w0^2/lambda; % rayleigh length
[XC,YC]=meshgrid(+(-R/2:R/2-1),-R/2:R/2-1);
theta=atan2(YC,XC);
ZoC=sqrt(XC.^2+YC.^2);
rou=ZoC*pixelsize/Magnify;
% z=[-1:0.1:1];% micron
zbar=z./z0;
wz=w0.*sqrt(1+zbar.^2);
U=zeros(R,R,length(z));
% U=(ramp(U,3))*1i+1;
for ii=1:length(z)
    roubar=rou./wz(ii);
    
    G=w0/wz(ii).*exp(-roubar.^2).*exp(1i.*roubar.^2.*zbar(ii)-1i*atan(zbar(ii)));
        [L]=Laguerre((n-abs(m))/2,abs(m),2*roubar.^2,R);
    RR=(sqrt(2).*roubar).^abs(m).*L;
%         theta=phiphi([R,R],'math');
        
    Phi=exp(1i*m.*theta);
    ZZ=exp(-1i*n*atan(zbar(ii)));
    
    Fig1=G.*RR.*Phi.*ZZ;
    U(:,:,ii)=Fig1./sqrt(sum(sum(Fig1.*conj(Fig1))));
    
end

