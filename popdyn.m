function [ dx ] = popdyn( t,x )
global r k be bp m u

% Apply penalty to carrying capacity due to evolved resistance
K=Kfun(u);
% Compute mu
mu=m./(k+be+bp*u);
% Compute population change
dx=x.*(r*((K-sum(x))./K)-mu);
end
