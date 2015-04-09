function [ dy ] = rabbitfox( t,y )
global a b c d
dy=zeros(2,1);
dy(1) = a*y(1) - b*y(2)*y(1);
dy(2) = c*y(1)*y(2) - d*y(2);
end

