function [ res ] = testode( t,x )
    global y n
    y = y + x;
    n = n + 1;
    res=cos(t)*(1+x);
end
