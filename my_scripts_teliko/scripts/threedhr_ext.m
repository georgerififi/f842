function f = threedhr_ext(X,I,a,b,c,d,r,s,xrest)

% Values of parameters
%a = 1; b = 3; c = 1; d = 5;
%r = .001;
%s = 4;
%xrest = -1.6;
%I = 5.5;

% fiducial point, i.e. the center of the sphere
x = X(1); 
y = X(2); 
z = X(3);

% orthonormal vectors of the ellipsoid
Y = [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];

f=zeros(12,1);

% 3dHR equation
%f(1) = y - a*x^3 + b*x^2 + (I-z);
%f(2) = c - d*x^2 - y;
%f(3) = r*(s*(x-xrest) - z);
f(1:3) = threedhr(X(1:3),I,a,b,c,d,r,s,xrest);

%Linearized system
Jac = [-3*a*x*2 + 2*b*x,  1, -1;
       -2*d*x,           -1,  0;
       r*s,               0, -r];
  
%Variational equation   
f(4:12) = Jac*Y;
