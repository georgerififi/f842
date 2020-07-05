function [lambda,t1,x1,x2] = lyapunov2(xinit,delta0,...
                    a,b,c,d,r,s,xrest,I,...
                    tstart,dt,T,...
                    plotflag,skiptransientflag)

%get time series of x y and z
%[x1,y1,z1] = Hindmarsh_Rose(P0,parameters,dt,T,tail,plotxyz);
tspan = tstart:dt:T;
[t1,x1] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit);
threshold = log2(sqrt(x1(:,1).^2+x1(:,2).^2+x1(:,3).^2));
%save final values in case needing to scan through parameter space later

%get slightly different trajectory time series
%[x2,y2,z2] = Hindmarsh_Rose(P0+delta0,parameters,dt,T,tail,plotxyz);
[t2,x2] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,xrest), tspan, xinit+delta0);
%save time series if opted for

%get squared distance
delta = ((x1(:,1)-x2(:,1)).^2 +(x1(:,2)-x2(:,2)).^2 +(x1(:,3)-x2(:,3)).^2);

%get plotting values of distance
delta = log2(sqrt(delta));

%find when delta reaches magnitude of attractor
i_leveled_out = find(delta>threshold);
if isempty(i_leveled_out)
  i_leveled_out = length(delta);
end
i_leveled_out = i_leveled_out(1);
%{
%determine sufficient time after transience
ten_percent_of_time = floor(0.1*i_leveled_out(1));
y_bot = mean(delta(1:(1+ten_percent_of_time)));
x_bot = mean(1:(1+ten_percent_of_time));
%get approximated value for maximum distance to consider
y_top = mean((delta((i_leveled_out(1)-ten_percent_of_time):(i_leveled_out(1)))));
x_top = mean((i_leveled_out(1)-ten_percent_of_time):(i_leveled_out(1)));
%calculate slope/ Lyapunov exponent
lambda = (y_top-y_bot)/(dt*(x_top-x_bot));
%}

ytofit = delta(1:i_leveled_out);
xtofit = t1(1:i_leveled_out);
%deltavalid = delta(valid);
%tvalid = t1(valid);
w = polyfit(xtofit, ytofit, 1);
lambda = w(1);
linearfit = w(1)*xtofit+ w(2);
%plot divergence time series if opted for
if plotflag==1
    figure;
    plot(xtofit,linearfit);
    hold on; plot(t1, delta);
    plot(t1,threshold,'k--');
    line([t1(i_leveled_out), t1(i_leveled_out)],...
        [-5, 5], ...
        'Color','black', ...
        'LineStyle', '--');
    %plot(dt*(1:(1/dt):T),delta(1:(1/dt):T));  hold on;  
    %plot(dt*[x_bot x_top],[y_bot y_top],'k');
    %plot(tvalid,deltavalid); hold on;
    %plot(tvalid,linearfit);
    %plot(tvalid,threshold(valid));
    xlabel('time'); ylabel('log2(|\delta|)');
    title(['\lambda = ',num2str(lambda)],'Interpreter','tex' );
end
end