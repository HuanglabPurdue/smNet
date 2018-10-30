function genfigs(obj)
zpos = obj.Zpos;
Num = numel(zpos);
psf = obj.ScaledPSFs;
xN = round(sqrt(Num));
yN = ceil(Num/xN);
h = figure;
h.Position = [200,300,100*xN,105*yN];
for n=1:Num
    if (n<=xN*yN)
        ha=axes;
        ii=mod(n,xN);
        jj=floor(n/xN)+1;
        if ii==0
            ii=xN;
            jj=n/xN;
        end
        psfi = psf(:,:,n);
        ha.Position=[(ii-1)/xN,(yN-jj)/yN,1/xN,1/yN];
        imagesc(psfi,[min(psfi(:)),max(psfi(:))]);axis off;axis equal;
        text(2,3, ['z=',num2str(zpos(n),3),'\mum'],'color',[1,1,1],'fontsize',12);
    end
end
colormap(jet);
h.InvertHardcopy = 'off';
h.Position = [200,300,105*xN,105*yN];
end
