clear; clc; close all;
a = 1; b = 3; c = 1; d = 5;
% the equilibrium points that come from the intersection 
% of nullclines
xlimt = 3;
x0 = -xlimt:.01:xlimt;

trace = -3*a*x0.^2 + 2*b*x0 - 1;
detrm = 3*a*x0.^2 + 2*(d-b)*x0;
indcs_t = findzcs(trace);
indcs_d = findzcs(detrm);

% figure of table 1
figure; 
p(1) = plot(x0,trace, 'DisplayName', 'Trace', ...
    'Color','blue', 'LineStyle', '--');  hold on;
p(2) = plot(x0,detrm, 'DisplayName', 'Determinant', ...
    'Color', 'red', 'LineStyle', '--'); 
p(3) = plot(x0,zeros(length(x0),1), 'Color', 'black');
for i = 1:length(indcs_t)
    line([x0(indcs_t(i)), x0(indcs_t(i))], ...
        [-10, 10], 'Color', 'blue', ...
        'DisplayName', sprintf('%.2f',x0(indcs_t(i)))); 
    hold on;
end
for i = 1:length(indcs_d)
    line([x0(indcs_d(i)), x0(indcs_d(i))], ...
        [-10, 10], 'Color', 'red', ...
        'DisplayName', sprintf('%.2f',x0(indcs_d(i))));
    hold on;
end
xlabel('x_0');
set(get(get(p(3),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%legend([p(1) p(2)]);
legend show
grid on; 
set(gca,'xtick',[-xlimt:.5:xlimt])
ylim([-10,10]);
title(sprintf('5 regions, a=%.2f, b=%.2f, c=%.2f, d=%.2f', ...
    a,b,c,d));

% figure of the equilibrium soltns: q = x^3 + px^2
p = (d-b)/a;
q = c/a;
fx = x0.^3 + p*x0.^2;
figure;
plot(x0, fx, 'DisplayName', 'x^3+px^2'); hold on;
line([-xlimt, xlimt], [q,q], ...
    'LineStyle', '--', 'DisplayName', 'q', 'Color', 'red');
line([-xlimt, xlimt], [(4*p^3)/27,(4*p^3)/27], ...
    'LineStyle', '--', 'DisplayName', '4p^3/27', 'Color', 'black');
grid on;
legend show;
xlabel('x_0')
title('equilibrium points x_0, q=x_0^3 + p x_0^2')

% check that they are the same as the zero crossings
D0L = -2*(d-b)/(3*a); % 0 cross det left
D0R = 0; % 0 cross det right
T0L = (b-sqrt(b^2-3*a))/(3*a);
T0R = (b+sqrt(b^2-3*a))/(3*a);

fprintf('DOL=%.2f, DOR=%.2f, TOL=%.2f, TOR=%.2f\n', ...
    D0L, D0R, T0L, T0R)

%helper
function indcs = findzcs(v) % find zero crossings
epsilon=.01;
v(v==0) = epsilon;
indcs = find(v(:).*circshift(v(:),-1) < 0);
end