function y = rk( t,x,h,treat )
    k1=h*dosedyn(t,x,treat(t));
    k2=h*dosedyn(t+h/2,x+k1/2,treat(t+h/2));
    k3=h*dosedyn(t+h/2,x+k2/2,treat(t+h/2));
    k4=h*dosedyn(t+h,x+k3,treat(t+h));
    y=x+(k1+2*(k2+k3)+k4)/6;
end
