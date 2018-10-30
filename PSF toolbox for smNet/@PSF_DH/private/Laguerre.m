% Laguerre function
function [L]=Laguerre(n,alpha,x,R)
L=zeros(R,R);
La=ones(R,R);
Lb=(-x+alpha+1).*ones(R,R);
if n==0, L=La;
else if n==1, L=Lb;
    else
        for ii=2:n
            L=(2*ii-1+alpha-x)/ii.*Lb-(ii-1+alpha)/ii.*La;
            La=Lb;
            Lb=L;
        end
    end
end

