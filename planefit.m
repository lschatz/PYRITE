function f = planefit(x,y,z)

%f = planefit(x,y,z)
%
% Fit a plane to the function defined by vectors x, y, z
% Fit is given by z = f(1)*x+f(2)*y+f(3)

m = zeros(3);
v = zeros(3,1);
m(1,1) = sum(x.*x);
m(1,2) = sum(x.*y);
m(1,3) = sum(x);
m(2,1) = m(1,2);
m(2,2) = sum(y.*y);
m(2,3) = sum(y);
m(3,1) = m(1,3);
m(3,2) = m(2,3);
m(3,3) = length(z);
v(1) = sum(z.*x);
v(2) = sum(z.*y);
v(3) = sum(z);
inv = invert3x3(m);
f = inv*v;

