function lam = lyapunov3(x,dt)
% https://www.math.tamu.edu/~mpilant/math442/Matlab/examples.html
% Wolf, 1985 section 5.1
% calculate lyapunov coefficient of time series

ndata=size(x,1);

%N2 = floor(ndata/2);
%N4 = floor(ndata/4);
N = ndata;
N3 = floor(ndata/3);
M = 10; % lenght of trajectory to study divergence
TOL_norm = 1.0e-8; % fiduciary and test point should have a norm-difference at least TOL_norm. 
TOL_time = 5; % fiduciary and test point should be at least TOL_time sample points apart.

exponent = nan(N-2*N3,1); %(N4+1,1);

for i=N3:(N-M)%(N-N3)  % middle third of the data should be 'not-transient'
   if mod(i,500)==0; fprintf('.'); end % visualize evolution
   dist = norm(x(i+1,:)-x(i,:));
   indx = i+1;
   for j=N3:(N-M)%((N-N3)-M) %1:(ndata-len_traj)
       if (i ~= j) && ...
            norm(x(i,:)-x(j,:))<dist && ...
            norm(x(i,:)-x(j,:))>TOL_norm && ... % norm-space
            ((j>i+TOL_time) || (j<i-TOL_time)) % time-space
                dist = norm(x(i,:)-x(j,:));
                indx = j; % closest point!
       end
   end
   % if, after all the searching, the found closest point is at leat 
   % TOL in norm-space and TOL_time in time away from point 'i', 
   % then estimate local rate of expansion (i.e. largest eigenvalue) for length M.
   %if norm(x(i,:)-x(indx,:))>TOL_norm && ... % norm-space
   %   ((indx>i+TOL_time) || (indx<i-TOL_time)) % time-space
   if indx~=i+1
       fprintf('point %d checkin , with %d \n',i, indx)    
       exponent(i-N3+1) = (log2(norm(x(i+M,:)-x(indx+M,:))) ...
                          - log2(norm(x(i,:)-x(indx,:)))) ...
                          / (M*dt); %  formula : (1/(tM-t0)) * sum_k={1..M} (L(i+k)/L(i+(k-1)))
   end
end

figure; plot(exponent); title('λ'); % plot the estimates for each initial point (fairly noisy)

% now, calculate the overal average over N4 data points ...
lam = nanmean(exponent); % return the average value
% if lam > 0, then system is chaotic
