function [ du ] = stratdyn( t,u )
global r bp m Kmax x sig
X=sum(x);
du=(r*X*exp(u.^2./(2*sig^2)).*u./(sig^2).*Kmax)-m*bp;
end
