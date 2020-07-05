% PARAMETERS
clear;  clc;  close all;
a = 1.0; b = 3; c = 1.0; d = 5.0;
r = .001; %.003;
s = 4; % 1; adaptation / 4; bursting
xrest = -1.6; % resting potential
xinit = [-1.618;-12.09;0]; % the left eqlbm (stable node)
dt = .05; % time step of numerical integration
T = 3000; % length for numerical integration
tspan = 0:dt:T; % time span
%tskip = ceil(T/3);
t_defaultskipidx = 3000; % default skip index
%eps = 1e-5;
Imin = 1;
Imax = 3.5; % max val of I to try
Istep_bif = .005; % I step for bif (fig.2)
Ilist = Imin:Istep_bif:Imax; % I vals for bifurcation (fig.2)
Istep_ts = .5; % I step for time series fig.1
Ilistts = Imin:Istep_ts:Imax;
n_subplots = length(Ilistts); % # time series subplots of fig.1
xthreshmax = 1.4; % xmax > 1.6: bursting peaks. i dont use it because i convolve to denoise, values change...
xthreshmin = -.9; % xmin > -1: bursting, xmin < -1: between bursts. i dont use it because i convolve to denoise, values change...
xthreshminrec = -1.3; % threshld for x min in recovery period 
%Tthresh =  120; % decide by observation of the recovery period duration (fig.1 : time series for different I vals)
kernel = [1 1 2 3 2 1 1];
kernel = kernel./sum(kernel); % kernel for denoising

% run exprmnt
k=1;
kk=1;
for I = Ilist  
    if (mod(I,.5)==0); fprintf('.'); end
    xmax = []; % burst peak vals
    xmin = []; % 'burst' lows and 'in-between-bursts' lows
    tstmps = []; % burst peaks timestamps
    
    % run 3D HR ode
    %f = @(t,x) [x(2) - a*x(1)^3 + b*x(1)^2 + (I - x(3)); ...
    %            c - d*x(1)^2 - x(2); ...
    %            r*(s*(x(1)-xrest) - x(3))];
    %[t,x] = ode45(f, tspan, xinit);
    [t,x] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit);
    
    %%%%%
    % fig.1. plot time series
    %%%%%
    if (mod(I,Istep_ts)==0)
        figure(1);
        subplot(n_subplots,1,k);
        plot(t,x(:,1));
        grid on;
        ylim([-2, 3]);
        ylabel('x')
        xlabel('t')
        title(sprintf('I:%.1f, s:%.1f',I,s));
        k=k+1;
    end
    
    %%%%%
    % fig.2 bifurcation diagram: plot (I,x), for x=xmin 
    % (either recovery or burst lows: black) and x=xmax (burst: red)
    %%%%%
    
    % denoise and keep only membrane potential
    x = conv(x(:,1),kernel,'same'); 
    
    % skip values
    indcs = find(x(t_defaultskipidx:end) < xthreshminrec); 
    if isempty(indcs); tskip = t(t_defaultskipidx); % default skip if doenst have recovery
    else; tskip = t(indcs(1)+t_defaultskipidx); % time stamp of first recovery
    end
    x = x(t>tskip); % skip everything BEFORE the first recovery 
    fprintf('skipped %.1f seconds \n',tskip)
    
    % update current hyperparameters
    n = length(x);
    %xthreshmax = prctile(x(x>0),95); % thresh to identify peaks. couldn't tune it so i turned back to hardocding it above.
    %xthreshmin = prctile(x(x<0),50); % thresh to identify lows. (...)
    diffx = diff(x);
    epsmax = prctile(diffx(diffx>0),5); % small epsilon. helps to identify peaks. suppose numerical noise 5%
    epsmin = prctile(diffx(diffx<0),95); % small epsilon. helps identify lows. suppose numerical noise 5%
    for i = 2:n-1 
        if (x(i) > (x(i-1)+epsmax)) && ... % maxs
            (x(i) > (x(i+1)+epsmax)) && ...
            (x(i) > xthreshmax)
            xmax = [xmax, x(i)]; % append
            tstmps = [tstmps, t(i)]; % append
        elseif (x(i) < (x(i-1)+epsmin)) && ... % mins
                (x(i) < (x(i+1)+epsmin)) && ...
                (x(i) < xthreshmin)
            xmin = [xmin, x(i)]; % append
        end 
    end
        
    figure(2)
    scatter(I*ones(length(xmax),1),xmax,'r.');
    scatter(I*ones(length(xmin),1),xmin,'k.');
    hold on; grid on;
    
    %%%%%
    % fig.3. plot period bifurcation diagram: (I,T), for 
    % T=Tburst_min, or Tburst_max (red) and T=Trecovery (black) 
    %%%%%
    
    Tbetweenpeaks = diff(tstmps);
    %Tbetweenpeaks_burst = Tbetweenpeaks(Tbetweenpeaks<Tthresh);
    %Tbetweenpeaks_inbetween =  Tbetweenpeaks(Tbetweenpeaks>=Tthresh);
    figure(3);
    scatter(I*ones(length(Tbetweenpeaks),1), ...
            Tbetweenpeaks,'b.');
    %{
    scatter(I*ones(length(Tbetweenpeaks_inbetween),1), ...
            Tbetweenpeaks_inbetween,'b.');
    Tmin = min(Tbetweenpeaks_burst);
    Tmax = max(Tbetweenpeaks_burst);
    if ~isempty(Tmin); scatter(I,Tmin(1),'r.'); end
    if ~isempty(Tmax); scatter(I,Tmax(1),'g.'); end
    %}
    hold on;
    grid on;
end 
figure(2); xlabel('I'); ylabel('x_{max},x_{min}'); title(sprintf('Bifurcation diagram for I, s:%.1f',s)); 
figure(3); xlabel('I'), ylabel('T');  set(gca,'yscale','log'); title(sprintf('Bifurcation diagram for I, s:%.1f',s));