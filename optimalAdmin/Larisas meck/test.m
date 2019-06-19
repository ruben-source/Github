% t = 5;
% y = 2*ones(1,3);
lambda = 0*ones(1,3);  obs_start = 1; obs_end = 30;
% adjfunc(t,y,lambda,g,obs_start, obs_end)
% Forwardfunc(y)
% ForwardfuncJacobi(y)
% A=AdjfuncJacobi(y)

time_mesh = 0:20;
alpha = 10^-9*ones(5,length(time_mesh));
x_initial = [1; 1; 1];
g = zeros(3,length(time_mesh));
x = ForwardNewton(alpha,time_mesh,x_initial);
lambda = AdjointNewton(alpha, x, g, time_mesh, obs_start, obs_end);
plot(time_mesh,x)



%%
gn = gradientfunc(s-1)
gm = gradientfunc(s)
sigma = norm(gm)/norm(gn);
d = -gm + sigma*d
















