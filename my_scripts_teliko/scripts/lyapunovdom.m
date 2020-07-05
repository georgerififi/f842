function [lambda,t1,x1,x2,delta,threshold] = lyapunovdom(P0,delta0,...
                                            a,b,c,d,r,s,x0,I,...
                                            tstart,tstep,tend)                

tspan = tstart:tstep:tend;

[t1_wt,x1_wt] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,x0), tspan, P0); % x1 with transient
[~,x2_wt] = ode45(@(t,x) threedhr(x,I,a,b,c,d,r,s,x0), tspan, P0+delta0); % x2 with transient
% assert(t1==t2)

% 20% transient transient: throw it away
twopc = ceil(length(tspan)/5); 
t1 = t1_wt(twopc:end);
clear t1_wt
x1 = x1_wt(twopc:end,:);
clear x1_wt
x2 = x2_wt(twopc:end,:);
clear x2_wt

threshold = log2(sqrt(x1(:,1).^2 + x1(:,2).^2 + x1(:,3).^2));
%threshold2 = log2(sqrt(x2(:,1).^2 + x2(:,2).^2 + x2(:,3).^2));
%threshold = (threshold1+threshold2)/2;

%get squared distance
delta = ((x1(:,1)-x2(:,1)).^2 ...
        +(x1(:,2)-x2(:,2)).^2 ...
        +(x1(:,3)-x2(:,3)).^2);

%get plotting values of distance
delta = log2(sqrt(delta));

% smooth
discret = 500;
len = ceil((tend/discret) / tstep);
kernel = (1/len)*ones(1,len);
delta_smooth = conv(delta, kernel, 'same');
threshold_smooth = conv(threshold, kernel, 'same');

%find when delta reaches magnitude of attractor
i_leveled_out = find(delta_smooth > threshold_smooth);
if isempty(i_leveled_out)
  i_leveled_out = length(delta);
end
i_leveled_out = i_leveled_out(1);


%determine sufficient time after transience
ten_percent_of_time = floor(0.1*i_leveled_out(1));
y_bot = mean(delta(1:(1+ten_percent_of_time)));
x_bot = mean(1:(1+ten_percent_of_time));

%get approximated value for maximum distance to consider
y_top = mean((delta((i_leveled_out(1)-ten_percent_of_time):(i_leveled_out(1)))));
x_top = mean((i_leveled_out(1)-ten_percent_of_time):(i_leveled_out(1)));

%calculate slope/ Lyapunov exponent
lambda = (y_top-y_bot)/(tstep*(x_top-x_bot));

%ytofit = delta(1:i_leveled_out);
%xtofit = t1(1:i_leveled_out);
%w = polyfit(xtofit, ytofit, 1);
%lambda = w(1);
%linearfit = w(1)*xtofit+ w(2);


%plot divergence time series if opted for
subplot(211);
p(1) = plot(tstep*[x_bot+twopc x_top+twopc], [y_bot y_top], 'k'); hold on;
p(2) = plot(t1, delta);
p(3) = plot(t1,threshold);
%plot(t1, deltasmooth, 'r--')
%plot(dt*(1:(1/dt):T),delta(1:(1/dt):T,'r')); 
%plot(dt*[x_bot x_top],[y_bot y_top],'k');
%plot(tvalid,deltavalid); hold on;
%plot(tvalid,linearfit);
%plot(tvalid,threshold(valid));
xlabel('time'); ylabel('log2(||\delta||)');
uistack(p(1),'top')

subplot(212); plot(t1,x1(:,1)); 
hold on; plot(t1,x2(:,1));
ylim([-2, 3])
%suplabel(sprintf('I: %.2f',I));