%%
clear;  clc;  close all;
a = 1.0; b = 3; c = 1.0; d = 5.0;
r = .001;
%s = 4; % 1; adaptation / 4; bursting
x0 = -1.6; % resting potential
xinit0 = [-6.5;0;0]; % below the separatrix -> converges to 
                     % the left eqlbm (stable node)
xinit1 = [-1.618;-12.09;0]; % the left eqlbm (stable node)
xinit2 = [-1.2;-2;0]; % above the separatrix -> converges 
                        % to the limit cycle 
xinit = [xinit0 xinit1 xinit2];
dt = .01;
T = 700;
t = 0:dt:T;
disp_step = ceil(length(t)/10);

%% (a) adaptation - fig.5. 
% In HIN84_2 we say how to keep the spiking going, even
% after the I step is off. This is biologically not possible.
% That's why we introduce the concept of adaptation, with the 
% 3d HR model.
s=1;
% define a I step (window, amplitude)
I = zeros(length(t),1);
tstamp = 30;
c_win = 70;
c_amp = 1;
I((tstamp<t)&(t<tstamp+c_win)) = c_amp;
figure; 
splot(1) = subplot(211); hold on; grid on;
for i = 1:size(xinit,2)
    x = nan(3,T);
    x(:,1) = xinit(:,i);
    for j=2:length(t)
        if mod(j,disp_step)==0; fprintf('.'); end
        x_new = rk(@threedhr,x(:,j-1),dt,I(j),a,b,c,d,r,s,x0);
        x(:,j) = x_new;
    end
    plot(t,x(1,:), ...
        'DisplayName',sprintf('x_{init}=(%.1f,%.1f,%.1f)',...
            xinit(1,i), xinit(2,i), xinit(3,i)));
end
line([t(1), t(end)],[x0, x0], ...
    'LineStyle', '--', 'Color', 'black', 'DisplayName', 'x_{rest}');
ylabel('x'); legend show;
splot(2) = subplot(212); 
plot(t,I,'Color', 'black'); 
ylabel('I'); xlabel('t')
linkaxes(splot, 'x')

%% (b) bursting - fig.6.
s=4;
xinit = [xinit1];
dt = .1;
T = 5000;
t = 0:dt:T;
disp_step = ceil(length(t)/10);
Ivals = [.4 2 4];
for ii = 1:length(Ivals)
    I = Ivals(ii);
    figure(ii+1);
    splot(1) = subplot(211); hold on; grid on;
    for i = 1:size(xinit,2)
        x = nan(3,T);
        x(:,1) = xinit(:,i);
        for j=2:length(t)
            if mod(j,disp_step)==0; fprintf('.'); end
            x_new = rk(@threedhr,x(:,j-1),dt,I,a,b,c,d,r,s,x0);
            x(:,j) = x_new;
        end
        plot(t,x(1,:), ...
            'DisplayName',sprintf('x_{init}=(%.2f,%.2f)',...
                xinit(1,i), xinit(2,i)));
    end
    ylabel('x'); legend show;
    splot(2) = subplot(212); 
    line([t(1) t(end)],[I I],'Color', 'black'); 
    ylabel('I'); xlabel('t'); 
    ylim([0,max(Ivals)]);
    linkaxes(splot, 'x')
end

%% (c) random structure
xinit = [xinit1];
dt = .1;
T = 5000;
t = 0:dt:T;
I = 3.25; % change I
r = .005; % change r
splot(1) = subplot(211); hold on; grid on;
for i = 1:size(xinit,2)
    x = nan(3,T);
    x(:,1) = xinit(:,i);
    for j=2:length(t)
        if mod(j,disp_step)==0; fprintf('.'); end
        x_new = rk(@threedhr,x(:,j-1),dt,I,a,b,c,d,r,s,x0);
        x(:,j) = x_new;
    end
    plot(t,x(1,:), ...
        'DisplayName',sprintf('x_{init}=(%.2f,%.2f)',...
            xinit(1,i), xinit(2,i)));
end
ylabel('x'); legend show;
splot(2) = subplot(212); 
line([t(1) t(end)],[I I],'Color', 'black'); 
ylabel('I'); xlabel('t'); 
ylim([0,max(Ivals)]);
linkaxes(splot, 'x')

%% (d) post inhibitory rebound - fig.8
xinit = [xinit1];
dt = .05;
T = 800;
t = 0:dt:T;
I = zeros(length(t),1);
tstamp = 50;
c_win = 200;
c_amp = -3;
I((tstamp<t)&(t<tstamp+c_win)) = c_amp;
figure;
splot(1) = subplot(211); hold on; grid on;
for i = 1:size(xinit,2)
    x = nan(3,T);
    x(:,1) = xinit(:,i);
    for j=2:length(t)
        if mod(j,disp_step)==0; fprintf('.'); end
        x_new = rk(@threedhr,x(:,j-1),dt,I(j),a,b,c,d,r,s,x0);
        x(:,j) = x_new;
    end
    plot(t,x(1,:), ...
        'DisplayName',sprintf('x_{init}=(%.2f,%.2f)',...
            xinit(1,i), xinit(2,i)));
end
ylabel('x'); legend show;
splot(2) = subplot(212); 
plot(t,I,'Color','black'); 
ylabel('I'); xlabel('t'); 
ylim([c_amp,0]);
linkaxes(splot, 'x')

