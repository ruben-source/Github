function f = optimizationAlgorithmOutline(time_mesh,alpha_guess)
    alpha(0) = alpha_guess;
    x0 = ForwardNewton(alpha(0));
    lambda0 = AdjointNewton(alpha(0));
    g(0) = gradientfunc(alpha(0),lambda(0),x(0));
    d(0) = -g(0);
    beta(0) = 1/gamma(0); % d(0) = - g(0)
    %beta = opt_step_length; % first step length is pre-set
    g(0) = sign(g0);
    alpha(1) = alpha(0) + beta(0)*g(0);
    for s = 1:max_iter-1
        x(s) = ForwardNewton(alpha(s));
        lambda(s) = AdjointNewton(alpha(s),x(s));
        g(s) = gradientfunc(alpha(s),x(s),lambda(s));
        sigma(s) = (norm(g(s))/norm(g(s-1)))^2;
        d(s) = -g(s) + sigma(s)*d(s-1);
        gamma(s) = gamma(0)/(s+1)^0.5;
        beta(s) = g(s)*d(s)'./(gamma(s)*d(s)*d(s)');
        if g(s) <= theta
            %fill in vacant spots 
            alpha(s+1:end) = alpha(s);
            disp('Has converged')
            break
        elseif norm(g(s)) >= norm(g(s-1))
            alpha(s+1:end) = alpha(s);
            disp('Grows abruptly')
            break
        elseif 0.999999*norm(alpha(s-1,:)) < norm(alpha(s,:)) && norm(alpha(s,:)) < 1.000001*norm(alpha(s-1,:))
            alpha(s+1:end) = alpha(s);
            disp('All alpha has converged')
            break
        else
            g(s) = sign(g(s));
            for t = nodes
               if alpha(s+1,t) <= interval_min
                   alpha(s+1,t) = interval_min;
               elseif alpha(s+1,t) >= interval_max
                   alpha(s+1,t) = interval_max;
               elseif b(s) < 0.01
                   alpha(s+1,t) = alpha(s,t) + beta(s,t)*g(s,t);
               else
                   if g(s) == -g(s-1)
                        beta_opt(j) = beta_opt(j)/2;
                   end
                   alpha(s+1,t) = alpha(s,t) - beta_opt(j)*g(s,t);
               end
            end
        end
    end
    
end