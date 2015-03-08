function [ dx ] = tumor( t,x )
global r sig k be bp m Kmax v

% Apply penalty to carrying capacity due to evolved resistance
K=Kmax*exp((-v^2)/(2*sig^2));
% Growth inhibition due to population size
pen=(K-x)/K;
% Compute mu
mu=m/(k+be+bp*v);
dx=r*x*pen-mu*x;

end
