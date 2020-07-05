function [lambda,t1,x1,x2] = lyapunov2(P0,delta0,...
                    a,b,c,d,r,s,xrest,I,...
                    tstart,dt,T,...
                    plotflag, skiptransientflag)
% get time series of x y and z
tspan = tstart:dt:T;
xinit = P0;
[t1,x1] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit);

% get slightly different trajectory time series
xinit = P0 + delta0;
[t2,x2] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit);

%{
% if you think the transient is too big and ends up not using 
% the attractor area, then use this (?)
if skiptransientflag==1
    idx = find(x1(:,1)<-1.4);
    tskip = t1(idx(1));
    x1 = x1(t1>tskip,:);
    t1 = t1(t1>tskip);
    x2 = interp1(t2,x2,t1);
end
%}

% assert t1=t2
%all(t1==t2) 

% get squared distance
delta = log10(sqrt( (x1(:,1)-x2(:,1)).^2 ...
                  + (x1(:,2)-x2(:,2)).^2 ...
                  + (x1(:,3)-x2(:,3)).^2) );

% find when delta reaches magnitude of attractor
threshold = log10(sqrt(x1(:,1).^2 + x1(:,2).^2 + x1(:,3).^2));
i_leveled_out = find(delta>threshold);
if isempty(i_leveled_out); i_leveled_out = length(delta);
else; i_leveled_out  = i_leveled_out(1); % keep the first time this occurs. we only can use this
end

%{
% determine sufficient time after transience
ten_percent_of_time = floor(0.1*i_leveled_out);
y_bot = mean(delta(1:(1+ten_percent_of_time)));
x_bot = mean(1:(1+ten_percent_of_time));

% get approximated value for maximum distance to consider
y_top = mean((delta((i_leveled_out-ten_percent_of_time):(i_leveled_out))));
x_top = mean((i_leveled_out-ten_percent_of_time):(i_leveled_out));

%calculate slope/ Lyapunov exponent
lambda = (y_top-y_bot)/(dt*(x_top-x_bot));
%}

valid = (delta(1:i_leveled_out)~=Inf);
deltavalid = delta(valid);
tvalid = t1(valid);
w = polyfit(tvalid, deltavalid, 1);
lambda = w(1);
linearfit = w(1)*tvalid + w(2);

%plot divergence time series if opted for
if plotflag==1
    figure;
    %plot(dt*(1:(1/dt):T),delta(1:(1/dt):T));  hold on;  
    %plot(dt*[x_bot x_top],[y_bot y_top],'k');
    plot(t1,delta); hold on;
    plot(tvalid,linearfit);
    plot(t1,threshold);
    xlabel('t'); ylabel('log10(|\delta|)'); 
    title(['\lambda = ',num2str(lambda)],'Interpreter','tex' );
end
end