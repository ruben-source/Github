function f = forwardfunc(t, y, par, L_t,L_d, M_t,M_d, ...
                                    I_t,I_d, time_mesh)


%dot T = a*T*(1 - b*T) - c*N*T - D*T - K_T*(1 - exp(-M)*T);
%dot N = e*C - f*N + g*T^2*N/(h+^T^2) - p*N*T;
%dot L = -m*L + j*D^2*T^2/(k+D^2*T^2)*L - q*L*T + (r_1*N + r_2*C)*T -
%           u*N*L^2 - K_L*(1 - exp(-M))*L + p_I*L*I/(g_I+I) + v_L(t)
%dot C = alpha - beta*C -K_C*(1-exp(-M))*C
%dot M = -gamma*M + v_M(t)
%dot I = mu_I*I + v_I(t);
% D = d*(L/T)^l/(s+(L/T)^l)
% y(1) = T; y(2) = N; y(3) = L; y(4) = C; y(5) = M; y(6) = I;

    a = par(1);
    b = par(2);
    c = par(3);
    d = par(4);
    e = par(5);
    l = par(6);
    f = par(7);
    g = par(8);
    h = par(9);
    j = par(10);
    k = par(11);
    m = par(12);
    q = par(13);
    p = par(14);
    s = par(15);
    r_1 = par(16);
    r_2 = par(17);
    u = par(18);
    K_T = par(19);
    K_N = par(20);
    K_L = par(21);
    K_C = par(22);
    alpha = par(23);
    beta = par(24);
    gamma = par(25);
    p_I = par(26);
    g_I = par(27);
    mu_I = par(28);

    
    
    v_L =@(t) drugAdministration(t,L_t, L_d,time_mesh);   
    v_M =@(t) drugAdministration(t,M_t, M_d,time_mesh);     
    v_I =@(t) drugAdministration(t,I_t, I_d,time_mesh); 
    
    T = y(1); N = y(2);L = y(3); C = y(4); M = y(5); I = y(6);    
    D = d*(L/T)^l/(s+(L/T)^l);
    
    f1 = a*T*(1 - b*T) - c*N*T - D*T - K_T*(1 - exp(-M))*T;

    f2 = e*C - f*N + g*T^2*N/(h + T^2) - p*N*T - ...
            K_N*(1-exp(-M))*N;

    f3 = -m*L + j*D^2*T^2/(k+D^2*T^2)*L - q*L*T + ... 
            (r_1*N + r_2*C)*T - u*N*L^2 - K_L*(1 - exp(-M))*L + ...
            p_I*L*I/(g_I+I) + v_L(t);

    f4 = alpha-beta*C - K_C*(1-exp(-M))*C;

    f5 = -gamma*M + v_M(t);

    f6 = -mu_I*I + v_I(t);
    
    f = [f1; f2; f3; f4; f5; f6];

end
