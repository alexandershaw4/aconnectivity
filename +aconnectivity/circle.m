function [x,y] = circle(n,radius)

if nargin < 2 || isempty(radius)
    radius=3;%just an example
end

theta=linspace(0,2*pi,n);%you can increase this if this isn't enough yet
x=radius*cos(theta);
y=radius*sin(theta);

