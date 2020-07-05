clear; clc; close all;
xinitmain = [-1.60435320914467; 
            -11.8748883483645; 
            -0.0179916610514969];
xinit0 = [.1;.1;.1];
xinit1 = [3; 3; 0];
xinit2 = [-15; -4; 0];
xinit3 = [-10;-15; 0];
xinitspp = [xinit0 xinit1 xinit2, xinit3];
Ilist = [0, ... % quiescence (1 stable node)
        0.75, ... % transient burst (1 limit cycle and 1 stable node)
        1.2, ... % transient burst (1 limit cycle and 1 stable spiral)
        1.269866, ... % regular burst (hopf1: two limit cycles)
        2, 3, ... % regular burst (two limit cycles)
        3.25, ... % irregular spiking (chaotic)
        5.3989871, ... % regular spiking (hopf2: 1 limit cycle)
        6, ... % regular spiking
        6.2028149]; % regular spiking (hopf3: ??)
a = 1.0; b = 3; c = 1.0; d = 5.0;
r = .001; s = 4;
xrest = -1.6; 
tstart = 0;
dt = .05;
tend = 10000;
tspan = tstart:dt:tend;

% run 3D HR ode's
for I = Ilist
    [t,x] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinitmain);
    figure; 
    subplot(121);
    plot3(x(:,1),x(:,2),x(:,3)); hold on; 
    scatter3(xinitmain(1),xinitmain(2),xinitmain(3),'r*');
    xlim([-2,3]); 
    zlim([0, 3]);
    xlabel('x'); ylabel('y'); zlabel('z'); 
    title(sprintf(['xyz evolution (I=%.2f), from:', vec_to_title(xinitmain)],I));
    subplot(122);
    i=1;
    for xinit = xinitspp
        [ttmp,xtmp] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit);
        p(i) = plot(xtmp(:,1),xtmp(:,2), 'DisplayName', sprintf('xinit: %.2f',xinit(1))); hold on;
        scatter(xtmp(end,1),xtmp(end,2), '*');        
        i=i+1;
    end
    title(sprintf('pplane, I=%.2f',I)); 
    xlabel('x'); ylabel('y')
    legend([p(1) p(2) p(3) p(4)]) %show;
    figure; plot(t,x(:,1)); ylabel('x'); xlabel('t'); 
    title(sprintf(['time series (I=%.2f), from:', vec_to_title(xinitmain)],I));
    ylim([-2,3]);
end