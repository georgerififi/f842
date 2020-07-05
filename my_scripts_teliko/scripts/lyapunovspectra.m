function [Texp,Lexp]=lyapunov(n,rhs_ext_fcn,fcn_integrator,tstart,stept,tend,ystart,ioutp);
% thanks to https://in.mathworks.com/matlabcentral/fileexchange/4628-calculation-lyapunov-exponents-for-ode
% n number of nonlinear odes, here Lorenz is 3D system, i.e. n=3. 

n2=n*(n+1); % total number of odes
nit = round((tend-tstart)/stept); % number of steps

% Memory allocation 
y = zeros(n2,1); 
cum = zeros(n,1);

% orthonormal frame of the initial condition 'ystart'
% x-axis: [1 0 0], 
% y-axis: [0 1 0], 
% z-axis: [0 0 1]
y(1:n) = ystart(:);
for i = 1:n
  y((n+1)*i)=1; 
end
t = tstart;

Lexp = [];
Texp = [];
% main loop
for ITERLYAP = 1:nit
  
  % Solution of extended ODE system from t to t+step, for initial condition 'y' 
  [~,Y] = feval(fcn_integrator,rhs_ext_fcn,[t t+stept],y);  

  t = t+stept; % update t for the next iteration
  y = Y(end,:); % get the LAST 12-dimensional value of Y, i.e. Y(t+tstep)

  % construct new orthonormal basis by gram-schmidt
  znorm = zeros(n,1);

  % v1' = v1 / |v1|, where v1 = [x(2) y(2) z(2)] i.o.w. [Y(4) Y(5) Y(6)]
  znorm(1) = sqrt(sum(y(n+1:2*n).^2));
  y(n+1:2*n) = y(n+1:2*n)/znorm(1);

  % % v2' = v2 - <v1',v2> v1' , where v2 = [x(3) y(3) z(3)] or [Y(7) Y(8) Y(9)]
  % % v2' = v3 - <v1',v3> v1' - <v2',v3> v2', where v3 = [x(4) y(4) z(4)] or [Y(10) Y(11) Y(12)]
  % ... 
  for i = 2:n
    for j = 1:(i-1)
      y(n*i+1:(i+1)*n) = y(n*i+1:(i+1)*n) - y(n*i+1:(i+1)*n) * y(n*j+1:(j+1)*n)' * y(n*j+1:(j+1)*n);  
    end
    znorm(i) = sqrt(sum(y(n*i+1:(i+1)*n).^2));
    y(n*i+1:(i+1)*n) = y(n*i+1:(i+1)*n)/znorm(i);
  end

  % update running vector magnitudes: cumulative and lp
  cum = cum + log2(znorm); % calculate
  lp = cum/(t-tstart); % normalize: Lexp = lim_{k->inf} {  (1/(k*dt)) * sum_{k:0..inf} log2(norm1, norm2, norm3)  } 

  % Keep track of lp values and t values.
  Lexp = [Lexp; lp'];
  Texp = [Texp; t];

  % display result
  if mod(ITERLYAP,ioutp)==0
    fprintf('t=%6.4f',t);
    for k=1:n 
      fprintf(' %10.6f',lp(k)); 
    end
    fprintf('\n');
  end

end

% plot the evolution
p = plot(Texp,Lexp);
xlabel('Time'); 
title(sprintf('Lyapunov exponents, with final values \\lambda_1=%.3f, \\lambda_2=%.3f, \\lambda_3=%.3f', ...
												Lexp(end,1), Lexp(end,2), Lexp(end,3)));
legend(p, sprintf('\\lambda_1'), sprintf('\\lambda_2'), sprintf('\\lambda_3'));