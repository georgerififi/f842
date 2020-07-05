function lyapprox(init1,init2,T)
% find the approximate value of the maximal lyapunov
% exponent by integrating ode forwards in time
% with two different (but close) coefficients
% over the time interval [0,T] using 100 time steps ...

% decrease relative tolerance for step size (default = 1.0e-3
options=odeset('RelTol',1.0e-7);

% take 100 time steps;
dt = T/100;

% integrate lorentz equations with different initial conditions
[t,x]=ode45(@lorz,[0:dt:T],init1,options);
[t,y]=ode45(@lorz,[0:dt:T],init2,options);

% non chaotic equations (e.g. linear spring) 
% give very small exponent (1.0e-12)
% [t,x]=ode45(@myspring,[0:dt:T],init1,options);
% [t,y]=ode45(@myspring,[0:dt:T],init2,options);

% calculate array size
N=length(t);

N4=round(N/4); % first quartile
N2=round(N/2); % second quartile (median)

% calculate L2-norm of difference
for i=1:N
    z(i)=norm(x(i,:)-y(i,:)); % find L2 norm
    % z(i) = max(abs(x(i,:)-y(i,:))); % find max norm
end

% normalize differences
z=z'/z(1);

% restrict interval to second quartile of data
a=t(N4:N2);
b=log2(z(N4:N2));

% plot log(|d(t)|/|d(0)|) vs t over range.
plot(a,b)

% the slope of this curve is the maximal lyapunov exponent
res=polyfit(a,b,1);
lyapunov_exponent=res(1)

function xdot = myspring(t,x)
% right hand side of 2x2 system of equations
% which govern the spring mass system
% x''+omega^2 x = 0

omega=1;
xdot = zeros(2,1);
xdot(1) = x(2);
xdot(2) = -omega^2 * x(1);

function xdot = lorz(t,x)  
  xdot = zeros(3,1);
  sig = 10.0;
  rho = 28.0;
  bet = 8.0/3.0;
  xdot(1) = sig*(x(2)-x(1));
  xdot(2) = rho*x(1)-x(2)-x(1)*x(3);
  xdot(3) = x(1)*x(2)-bet*x(3);


