function [ K ] = Kfun( v ) % Compute carrying capacity based on strategy
global Kmax sig
K=Kmax.*exp((-v.^2)./(2*sig^2));
end
