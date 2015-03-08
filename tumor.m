function [ dx ] = tumor( t,x )
global r sig k be bp m Kmax v
% Compute K
K=Kmax*exp((-v^2)/(2*sig^2));
% Compute mu
pen=(K-x)/K;
dx=r*x*pen-mu*x;

end
