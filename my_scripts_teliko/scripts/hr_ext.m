function f=hr_ext(t,X,a,b,c,d,r,s,x0,I)
% Copyright (C) 2004, Govorukhin V.N.

% Values of parameters
%a=1; b=3; c=1; d=5; x0=-1.6; I=3.25;
%r=.001; s=4;

% fiducial point, i.e. the center of the sphere
x=X(1); y=X(2); z=X(3);

% orthonormal vectors of the ellipsoid
Y= [X(4), X(7), X(10);
    X(5), X(8), X(11);
    X(6), X(9), X(12)];

f=zeros(12,1);
%Lorenz equation
f(1) = y - a*x^3 + b*x^2 +I-z;
f(2) = c - d*x^2 -y;
f(3) = r*(s*(x-x0)-z);

%Linearized system
Jac = [-3*a*x^2 + 2*b*x  1  -1;
        -2*d*x           -1  0;
        r*s              0   -r];
  
%Variational equation   
f(4:12) = Jac*Y;
