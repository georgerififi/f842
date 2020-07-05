function dxdt = twodhr(x,I,a,b,c,d)
dxdt = [x(2) - a*x(1)^3 + b*x(1)^2 + I; ...
           c - d*x(1)^2 - x(2)];
end
