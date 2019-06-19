function f = jacobian(y,par)
    T = y(1); N = y(2) ; L = y(3); C = y(4); M = y(5); I = y(6);
    kvot = (L/T)^l;
    Dden = s + kvot; %Denominator of D
    
    D = d*kvot/Dden;
    dDdT = -D*l/T + D^2/d*l/T;
    dDdL = D*l/L - D^2/d*l/L;
    %d(D^2T^2/(k+D^2T^2))dT = dD2dT
    D2den = k + D^2*T^2; 
    dD2dT = (2*D*dDdT*T^2 + D^2*2*T)/D2den - D^2*T^2*(2*D*dDdT*T^2 + D^2*2*T)/D2den^2;
    %d(D^2T^2L/(k+D^2T^2))dT = dD2dL  
    dD2dL = (2*D*dDdL*T^2/Dden + D^2*T^2/D2den^2*2*D*dDdL*T^2)*L + D^2*T^2/D2den;
    
    [a*(1-b*T)-a*T*b-c*N-D-K_T*(1-exp(-M))
        
    ]
    
end