function [sse,model]=mainlobeFit(x,data)
R=size(data,1);
x0=x(1);
y0=x(2);
sigmax=x(3);
A=x(4);
[X,Y]=meshgrid(-(-R/2:R/2-1),-R/2:R/2-1);
model=A.*(exp(-(X-x0).^2./sigmax./sigmax).*exp(-(Y-y0).^2./sigmax./sigmax)+...
    exp(-(X+x0).^2./sigmax./sigmax).*exp(-(Y+y0).^2./sigmax./sigmax));
sse=sqrt(sum(sum((model-data).^2)));