%% Observations g, given (or made) and unchangeable
%made
m = 100; %number of points
t_final = 100; 
time_mesh = linspace(0,t_final,m);
x_initial = [1; 1; 1];
alpha_exact = [10^(-9); 10^(-10); 10^(-8); 10^(-10); 10^(-10)].*ones(5,m);
% use function rand to add noise
g_made = ForwardNewton(alpha0, time_mesh, x_initial) + rand(1,m);
AdjointNewton(alpha, x, g, time_mesh, obs_start, obs_end)
gradientfunc(x,lambda,alpha,gamma,s,n)

%%
%given observations = g and when they are observed (obs_start, obs_end)
%initial guess of alpha alpha_guess(t) 5 x length(time_mesh) matrix
%initial values of x = [x_T x_M1 x_M2]
function Opt = OptimizationAlgorithm(g, time_mesh,obs_start,...
                                    obs_end,alpha_guess, x_initial)
    Max_iter = 1000;
    num_of_alpha = size(alpha_guess,1);
    nodes = length(time_mesh);
    
    gamma0 = 0.1*ones(num_of_alpha,1);
    gamma =@(s,i) gamma0(i)/s^0.5;
    
    alpha       = zeros(Max_iter,num_of_alpha,nodes);
    beta        = zeros(Max_iter,num_of_alpha);
    beta_opt    = 0.0001*ones(num_of_alpha,nodes);
    grad        = zeros(Max_iter,num_of_alpha,nodes);
    sigma       = zeros(Max_iter,num_of_alpha);
    d           = zeros(Max_iter,num_of_alpha,nodes);
    
    alpha(1,:,:) = alpha_guess;
    x0 = ForwardNewton(alpha(1,:,:), time_mesh, x_initial);
    lambda0 = AdjointNewton(alpha(1,:,:),x0,g,time_mesh,obs_start,obs_end);
    for i = 1:num_of_alpha
       grad(1,i,:) = gradientfunc(x0,lambda0,alpha,1,i);
       d(1,i,:) = -grad(1,i,:);
       beta(1,i) = 1 / gamma(1,i);
       dir = sign(grad(1,i,:));
       alpha(2,i,:) = alpha(1,i,:) + beta_opt(i,:).*dir;
    end
    
    flag = ones(num_of_alphas,1); % if flag(i) = 1 then continue calculate grad for alpha(i)
    for s = 2:Max_iter
       xs = ForwardNewton(alpha(s,:,:),time_mesh,x_initial);
       lambdas = AdjointNewton(alpha(s,:,:), xs, g, time_mesh, obs_start, obs_end);
       for i = 1:num_of_alpha
           if flag(i)
               grad(s,i,:) = gradientfunc(xs, lambdas,alpha,s,i);
               sigma(s,i) = (norm(grad(s,i,:))/norm(grad(s-1,i,:)))^2;
               d(s,i,:) = -g(s,i,:) + sigma(s,i)*d(s-1,i,:);
               beta(s,i) = g(s,i,:)*d(s,i,:)'./(gamma(s,i)*d(s,i,:)*d(s,i,:)');

           end
       end
        
    end
end