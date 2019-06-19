function lambda = adjfunc(t,alpha,x,lambda,g,obs_start, obs_end)
    % alpha = (d_m1, d_m2, a_t1, a_t2, k_12)
     r = 0.93;
     beta_T = 3*10^9;
     d_m1 = alpha(1);
     d_m2 = alpha(2);
     a_t1 = alpha(3);
     a_t2 = alpha(4);
     beta_M = 9*10^8;
     sigma_m1 = 0.173;
     sigma_m2 = 0.173;
     k_12 = alpha(5);
     
     x_T = x(1); x_M1 = x(2); x_M2 = x(3);
     
     lambda1 = -lambda(1)*r*(1-2*x_T/beta_T) + lambda(1)*d_m1*x_M1 - lambda(1)*d_m2*x_M2...
               -lambda(2)*a_t1*x_M1*(1-(x_M1+x_M2)/beta_M) + lambda(2)*k_12*x_M1...
               -lambda(3)*a_t2*x_M2*(1-(x_M1+x_M2)/beta_M) - lambda(3)*k_12*x_M1...
               +(x_T-g(1))*z(t,obs_start, obs_end);
     lambda2 = -lambda(2)*a_t1*x_T*(1-(2*x_M1+x_M2)/beta_M) + lambda(2)*sigma_m1...
               +lambda(3)*a_t2*x_T*x_M2/beta_M - lambda(3)*k_12*x_T + lambda(1)*d_m1*x_T...
               +lambda(2)*k_12*x_T + (x_M1-g(2))*z(t,obs_start, obs_end);
     lambda3 = lambda(3)*sigma_m2 - lambda(3)*a_t2*x_T*(1-(x_M1+2*x_M2)/beta_M)...
               -lambda(1)*d_m2*x_T + lambda(2)*a_t1*x_T*x_M1/beta_M...
               +(x_M2-g(3))*z(t,obs_start,obs_end);
     lambda = [lambda1; lambda2; lambda3];
     
    function zeta = z(t,obs_start, obs_end)
       if t >= obs_start & t <= obs_end
           zeta = 1;
       else
           zeta = 0;
       end
    end
end