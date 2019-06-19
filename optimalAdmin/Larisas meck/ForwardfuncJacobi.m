function f = ForwardfuncJacobi(x,alpha)
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
    
    f11 = r*(1-2*x_T/beta_T)-d_m1*x_M1+d_m2*x_M2;
    f12 = -d_m1*x_T;
    f13 = d_m2*x_T;
    f21 = a_t1*x_M1*(1-(x_M1+x_M2)/beta_M)-k_12*x_M1;
    f22 = a_t1*x_T*(1-(2*x_M1+x_M2)/beta_M)-sigma_m1-k_12*x_T;
    f23 = -a_t1*x_T*x_M1/beta_M;
    f31 = a_t2*x_M2*(1-(x_M1+x_M2)/beta_M)+k_12*x_M1;
    f32 = -a_t2*x_T*x_M2/beta_M+k_12*x_T;
    f33 = a_t2*x_T*(1-(x_M1+2*x_M2)/beta_M)-sigma_m2;
    
    f = [f11 f12 f13;
         f21 f22 f23;
         f31 f32 f33];
    
end