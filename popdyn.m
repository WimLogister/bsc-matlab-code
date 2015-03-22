function [ dx ] = popdyn( t,x )
global r k be bp m sig Kmax
dx=zeros(2,1);
% Apply penalty to carrying capacity due to evolved resistance
K=Kfun(x(2,:));

X=sum(x(1,:));
% Compute mu
mu=m./(k+be+bp*x(2,:));
% Compute population change
dx(1)=x(1,:).*(r*((K-X)./K)-mu);
dx(2)=((r*X*exp(x(2,:).^2./(2*sig^2)).*x(2,:))./(sig^2).*Kmax)-m*bp;
end
