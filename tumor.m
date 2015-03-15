function [ dx ] = tumor( t,x )
global r k be bp m u

% Apply penalty to carrying capacity due to evolved resistance
K=Kfun(u);
% Growth inhibition due to population size
pen=(K-x)/K;
% Compute mu
mu=m/(k+be+bp*u);
dx=r*x*pen-mu*x;
end
