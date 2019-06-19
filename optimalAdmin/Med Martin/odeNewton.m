
%%
% t_final = 120;
% time_mesh = linspace(0,t_final,1000);
% L_t = [7];      L_d = [8*10^8];     
% I_t = [];       I_d = [];      
% M_t=0:14:90;    M_d = ones(1,length(M_t));
% par = mouse();
% y0 = [10^6; 5*10^4; 100; 1.1*10^7; 0; 0];
% [t,y] = odeNewton(par, L_t,L_d, M_t,M_d, ...
%                     I_t,I_d, time_mesh,y0);
% disp('done')
% plott(t,y,t_final);

function [time_mesh y] = odeNewton(par, L_t,L_d, M_t,M_d, ...
                                    I_t,I_d, time_mesh,y0)
    l = length(time_mesh);
    dt = time_mesh(2:end)-time_mesh(1:end-1);
    y = zeros(l,6);
    y(1,:) = y0;
    tol = 0.001;
    for i = 1:l-1
        w = y(i,:);    
        while abs(dw) < tol
            g = w - (y(i,:) + dt(i)*forwardfunc(time_mesh(i),y(i,:),...
                            par,L_t,L_d, M_t,M_d, I_t, I_d, time_mesh)');
            g_prim = eye(size(y0)) - dt(i)*jacobi();
            dw = 
        end
    end
end

% g(x) = 0 => Newton : x(k+1) = x(k) - g(x(k))/g'(x(k))