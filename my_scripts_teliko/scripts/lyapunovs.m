clear; clc; close all;
x_init = [.1; .1; .1];

%% Lexp, with Wolf orthonormalization algorithm

% display the current lyapunov values every 'iter' iterations
iter_disp = 10; 
% tspan
tstart = 0;
tend = 10000;
tstep = 0.5;
% parameters
a=1; b=3; c=1; d=5; x0=-1.6; 
r=.001; s=4;
I=5; % CHOOSE I
% run
[~,Lexp] = lyapunovspectra(3, @(t,x) hr_ext(t,x,a,b,c,d,r,s,x0,I), @ode45, ...
                    tstart, tstep, tend, x_init, iter_disp);

%% dominant lambda as divergence of two 
% predetermined close points fed as input.
tstep = .01;
tend=5000;
delta0 = [.00001; .00001; .000001];
Ilist = [0, .75, 1, ...
        1.099, 1.1, 1.11, ...
        1.2, 1.269866, 1.4, ...
        2, 2.5, ...
        3, 3.15, 3.25, 3.35, ...
        4.5];
i = 1;
for I = Ilist
    figure(i); 
    [lambda(i),t1,x1,x2,delta,thresh] = lyapunovdom(x_init,delta0,...
                                            a,b,c,d,r,s,x0,I,...
                                            tstart,tstep,tend);

    %figure; xlabel('x'); ylabel('y'); zlabel('z');
    %plot3(x1(:,1),x1(:,2),x1(:,3)); hold on; 
    %plot3(x2(:,1),x2(:,2),x2(:,3)); 
    %scatter3(x_init(1),x_init(2),x_init(3),'r*');
    %scatter3(x_init(1)+delta0(1),...
    %        x_init(2)+delta0(2),...
    %        x_init(3)+delta0(3),'r*');
    i = i+1;
end

% plot (lambdas,I)
lambda_thresh = 0.005;
indcs = lambda <= lambda_thresh;
figure; stem(Ilist(indcs), lambda(indcs)); hold on; 
stem(Ilist(~indcs), lambda(~indcs), 'r'); 
ylabel('lambda'), xlabel('I'); xticks(Ilist);

%% dominant lambda with fixed evolution time Wolf - FAILS
% https://www.math.tamu.edu/~mpilant/math442/Matlab/examples.html
tstep = .05;
tend = 2000;
I = 2;
tspan = tstart:tstep:tend;
[t,x] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit);
domlambda = lyapunov3(x,tstep)
figure; plot(t,x);

% TODOS:
% lyapunov by Wolf : https://www.mathworks.com/matlabcentral/fileexchange/48084-wolf-lyapunov-exponent-estimation-from-a-time-series
% http://sprott.physics.wisc.edu/chaos/lyapexp.htm
