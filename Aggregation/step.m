function [y]=krok(t,x,h)
% Runge-kutta method to discretize the ODE system (you can use this literally, except of that your expressions won't be vectors, but scalars)
k1=h*fde(t,x);
k2=h*fde(t+h/2,x+k1/2);
k3=h*fde(t+h/2,x+k2/2);
k4=h*fde(t+h,x+k3);
y=x+(k1+2*(k2+k3)+k4)/6;

