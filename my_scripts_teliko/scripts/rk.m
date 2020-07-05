function x_next = rk(f,x,dt,I,a,b,c,d,r,s,x0) % runge-kupta-4
if nargin<9 % 2d HR
	k1 = f(x,I,a,b,c,d);
	k2 = f(x+dt/2*k1,I,a,b,c,d);
	k3 = f(x+dt/2*k2,I,a,b,c,d);
	k4 = f(x+dt*k3,I,a,b,c,d);
	x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
else % 3d HR
    k1 = f(x,I,a,b,c,d,r,s,x0);
	k2 = f(x+dt/2*k1,I,a,b,c,d,r,s,x0);
	k3 = f(x+dt/2*k2,I,a,b,c,d,r,s,x0);
	k4 = f(x+dt*k3,I,a,b,c,d,r,s,x0);
	x_next = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end