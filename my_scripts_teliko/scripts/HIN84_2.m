%%
clear;  clc;  close all;
a = 1.0; b = 3; c = 1.0; d = 5.0; 
I = 0;
xinit0 = [-6.5; 0]; % converges to left (stable node) eq.: -1.6
xinit1 = [-1.618;-12.09]; % converges to left (stable node) eq.: -1.6
xinit2 = [-1.2;-2]; % converges to the limit cycle 
xinit = [xinit0 xinit1 xinit2];
dt = .005;
T = 400;
t = 0:dt:T;
disp_step = ceil(length(t)/10);

%%
figure; 
for i = 1:size(xinit,2)
    x = nan(2,T);
    x(:,1) = xinit(:,i);
    for j=2:length(t)
      x_new = rk(@twodhr, x(:,j-1),dt,I,a,b,c,d);
      x(:,j) = x_new;
    end
    subplot(311); hold on; grid on;
    plot(t,x(1,:), ...
        'DisplayName',sprintf('x_{init}=(%.2f,%.2f)',...
            xinit(1,i), xinit(2,i)));
    ylabel('x'); legend show;
    subplot(312); 
    plot(t,x(2,:)); grid on; hold on; ylabel('y')
end
subplot(313); 
line([t(1), t(end)],[I I], 'Color', 'black'); 
ylabel('I'); xlabel('t')


%% 
% this won't work because, in order for the oscillation to stop 
% the current step needs to stop, when the phase point is BELOW
% the saddle point!
I = zeros(length(t),1);
tstamp = 30;
c_win = 20;
c_amp = 1;
I((tstamp<t)&(t<tstamp+c_win)) = c_amp;
figure; 
subplot(211); hold on; grid on;
for i = 1:size(xinit,2)
    x = nan(2,T);
    x(:,1) = xinit(:,i);
    for j=2:length(t)
        if mod(j,disp_step)==0; fprintf('.'); end
        x_new = rk(@twodhr,x(:,j-1),dt,I(j),a,b,c,d);
        x(:,j) = x_new;
    end
    plot(t,x(1,:), ...
        'DisplayName',sprintf('x_{init}=(%.2f,%.2f)',...
            xinit(1,i), xinit(2,i)));
end
ylabel('x'); legend show;
subplot(212); 
plot(t,I,'Color', 'black'); 
ylabel('I'); xlabel('t')

%% ... so we should reset I, on the fly, in a very specific time window
% that will allow the phase point x, to stay on a limit cycle and not 
% go back to the stable node (the game is rigged)!
% find the line y=ax+b, that, if the phase point x_new is above, then
% it's ok to stop the I step. why is it ok? Because the separatrix of the saddle 
% won't be able to send x back to the stable node (-1-sqrt(5),-13-5*sqrt(5)). 
% Instead the phase point will be 'helped' towards the new limit cycle.
% the y=ax+b line is defined by the eigenvector-basin of the saddle point
% (-1,-4). 
v = [1,-.9]; % eigenvector
pt = [-1,-4]; % saddle point
atmp = v(2)/v(1); % a of the line y = ax+b
btmp = pt(2) - a*pt(1); % b of the line

% init current
I = zeros(length(t),1);
tstamp = 30; % time stamp to turn I 'on'
I(t>tstamp) = 1;
once_flag = 1;

% init x 
x = nan(2,T);
x(:,1) = xinit0; % study the effect of the I step on case xinit2 that didn't fire without the I step.

for j=2:length(t)
    if mod(j,disp_step)==0; fprintf('.'); end
    x_new = rk(@twodhr,x(:,j-1),dt,I(j),a,b,c,d);
    x(:,j) = x_new;
    if (once_flag==1) && (t(j)>tstamp) && (x_new(2)-atmp*x_new(1)>btmp) 
        %x_new
        I(j+1:end) = 0;
        once_flag = 0;
    end
end
figure;
s(1) = subplot(311);
plot(t,x(1,:), 'DisplayName',sprintf('x_{init}=(%.2f,%.2f)',...
        xinit(1,i), xinit(2,i)));
ylabel('x'); grid on;
s(2) = subplot(312); 
plot(t,x(2,:)); grid on;
ylabel('y')
s(3) = subplot(313); 
plot(t,I,'Color', 'black'); 
ylabel('I'); xlabel('t'); grid on;
linkaxes(s, 'x')