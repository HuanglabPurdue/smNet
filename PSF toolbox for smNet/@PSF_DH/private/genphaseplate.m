function [TR_phase,TR_mag,TR_phase0,TR_complex,k_z,k_r,Phi] = genphaseplate(R,n0,lambda,NA,pixelsize,Magnify,zRange,plotflag)
pixelsize0 = pixelsize/2;
z0=0;%micron
U1_1=GLmode(1,1,R,n0,lambda,NA,pixelsize0,Magnify,z0);
U3_5=GLmode(3,5,R,n0,lambda,NA,pixelsize0,Magnify,z0);
U5_9=GLmode(5,9,R,n0,lambda,NA,pixelsize0,Magnify,z0);
U7_13=GLmode(7,13,R,n0,lambda,NA,pixelsize0,Magnify,z0);
U9_17=GLmode(9,17,R,n0,lambda,NA,pixelsize0,Magnify,z0);
rotatebeam=U1_1+U3_5+U5_9+U7_13+U9_17;

% TR_complex=fftshift(fft2(rotatebeam)); % transfer function
TR_complex=double(ft(rotatebeam));
temp=abs(TR_complex);
temp(temp==0)=1e-20;
TR_phase0=TR_complex./temp;% phase plate

if plotflag == 1
    figure('position',[200,300,510,250],'color',[1,1,1])
    ht(1)=subplot('position',[0,0,1/2,1]);
    image(abs(TR_complex),'CDataMapping','scaled');
    ht(2)=subplot('position',[1/2,0,1/2,1]);
    image(angle(TR_complex),'CDataMapping','scaled');
    colormap(gray)
    axis(ht,'equal')
    axis(ht,'off')
    
end
%% generate GL mode around key mode

mx=[-11:1:11];
GLn=[];
GLm=[];
s=0;

for ii=1:length(mx)
    for jj=1:10
        GLn(jj+s)=abs(mx(ii))+2*(jj-1);
        GLm(jj+s)=mx(ii);
        if GLn(jj+s)>20
            break;
        end
    end
    s=s+jj;
end

keyGLm=[1 3 5 7 9];
keyGLn=[1 5 9 13 17];
cloudn=[];
cloudm=[];
g=1;
for ii=1:length(GLn)
    
    dmn=sqrt((GLn(ii)-keyGLn).^2+(GLm(ii)-keyGLm).^2);
    if min(dmn)<=2
       cloudn(g)=GLn(ii);
       cloudm(g)=GLm(ii);
       g=g+1;
    end

end 

if plotflag == 1
    
    figure('position',[200,300,500,400],'color',[1,1,1])
    set(gca,'fontsize',12)
    plot(GLm,GLn,'b.','Markersize',10)
    hold on
    plot(cloudm,cloudn,'go','linewidth',1.5)
    plot(keyGLm,keyGLn,'ro','linewidth',1.5)
    
    xlabel('m')
    ylabel('n')
    axis equal
    axis tight
    ylim([0,20])
end
%% generate GL mode in Fourier domain
Umn=zeros(R,R,length(cloudn));
FU=zeros(R,R,length(cloudn));
for k=1:length(cloudn)
   Umn(:,:,k)=GLmode(cloudm(k),cloudn(k),R,n0,lambda,NA,pixelsize0,Magnify,z0);
end

for k=1:length(cloudn)
    Fig1=squeeze(Umn(:,:,k));
    Fig3=double(ft(Fig1));
    FU(:,:,k)=Fig3./sqrt(sum(sum(conj(Fig3).*Fig3)));
end

%% mainlobe weight function generation
scale=R*pixelsize/Magnify;
Freq_max=NA/lambda;
[XC,YC]=meshgrid(-(-R/2:R/2-1),-R/2:R/2-1);
Phi=atan2(YC,XC);
ZoC=sqrt(XC.^2+YC.^2);
k_r=ZoC./scale;
NA_constrain=k_r<Freq_max;
k_z=sqrt((n0/lambda)^2-k_r.^2).*NA_constrain;
N=length(zRange);
mainlobe=zeros(R,R,N);
mainlobe1=zeros(R,R,N);
lobePos=[];
TR_mag=NA_constrain./sqrt(sum(sum(NA_constrain.^2))); %normalize
for j=1:N

    defocus_phase=2*pi*zRange(j).*k_z.*1i;
    pupil_complex=exp(defocus_phase).*TR_phase0.*TR_mag;
    Fig=fftshift(fft2(pupil_complex));
    Data=abs(Fig)./max(max(abs(Fig)));
    [max1,xi]=max(max(Data.^3,[],2));
    [max1,yi]=max(max(Data.^3,[],1));
    startpoint=[-yi+R/2,xi-R/2,2.5,1];
    estimate=fminsearch(@(x) mainlobeFit(x,Data.^2),startpoint,optimset('MaxIter',50,'Display','off'));    
    lobePos=cat(1,lobePos,estimate);   
    if zRange(j)==0
        a=j;
        break;
    end
end

for j=a-1:-1:1
lobePos=cat(1,lobePos,[lobePos(j,1),-lobePos(j,2),lobePos(j,3:4)]);
end

for j=1:N
    [sse, Fig1]=mainlobeFit(lobePos(j,:),Fig);
     mainlobe(:,:,j)=Fig1;% generate weight function
end

for j=1:N

    defocus_phase=2*pi*zRange(j).*k_z.*1i;
    pupil_complex=exp(defocus_phase).*TR_phase0.*NA_constrain;
    Fig=fftshift(fft2(pupil_complex));
    Fig1=maxf(abs(Fig),3)-medif(abs(Fig),8);
    mainlobe1(:,:,j)=Fig1;% generate weight function
 
end


%% iteration begin
TR_phase=TR_phase0;
for nn=1:10
    % initial phase plate decomposition
    
    TR_mag=NA_constrain./sqrt(sum(sum(NA_constrain.^2))); %normalize
    TR_complex=squeeze(TR_phase).*TR_mag;
    TRcoeff=[];
    
    for ii=1:length(cloudn)
        TRcoeff(ii)=real(sum(sum(TR_complex.*conj(squeeze(FU(:,:,ii))))));
    end
    
    % add rotation constrain
    Dmn=[];
    q=0.05;
    Mrotatebeam=zeros(R,R);
    MTR_phase=zeros(R,R);
    for k=1:length(cloudn)
        Dmn(k)=1;
        for g=1:length(keyGLn)
            Dmn(k)=((cloudm(k)-keyGLm(g))^2+(cloudn(k)-keyGLn(g))^2)^q*Dmn(k);
        end
    end
    MaxDmn=max(Dmn);
    
    for k=1:length(cloudn)
        Fig1=squeeze(Umn(:,:,k));
        Mrotatebeam=Fig1.*real(TRcoeff(k)).*(MaxDmn-Dmn(k))+Mrotatebeam;
    end
    
    for k=1:length(cloudn)
        Fig1=squeeze(FU(:,:,k));
        MTR_phase=Fig1.*real(TRcoeff(k)).*(MaxDmn-Dmn(k))+MTR_phase;
    end
    MTR_phase=MTR_phase./abs(MTR_phase);
    
    % add main lobe constrain
    psfA=zeros(R,R,N);
    TR_phasei=zeros(R,R,N);
    
    for j=1:N
        
        defocus_phase=2*pi*zRange(j).*k_z.*1i;
        pupil_complex=TR_mag.*exp(defocus_phase).*MTR_phase;
        Fig=fftshift(fft2(pupil_complex));
        psfA(:,:,j)=Fig;
        Fig1=squeeze(mainlobe(:,:,j));
        modpsfA=Fig.*Fig1;
        Fig2=fftshift(double(ft(modpsfA))).*exp(-defocus_phase);
        TR_phasei(:,:,j)=Fig2;
        
    end
    Fig3=squeeze(mean(TR_phasei,3));
    TR_phase=Fig3./abs(Fig3);
end
TR_phase = TR_phase.*NA_constrain;
end