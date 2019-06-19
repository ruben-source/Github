% T = y(:,1) = tumor cell population
% N = y(:,2) = total NK cell population,
% L = y(:,3) = total CD8+T cell population;
% C = y(:,4) = number of circulating lymphocytes (or white blood cells);
% M = y(:,5) =  chemotherapy drug concentration in the bloodstream
% I = y(:,6) =  immunotherapy drug concentration in the bloodstream
% v_L = TIL drugs given (Immunotherapy adding antigen-specific cytolytic immune cells)
% v_I =  amount of immunotherapy drug given 
% v_M =  amount of chemotherapy drug given 

%%
clc
disp('section 5.1: Immune system response to tumor: mouse data')
disp('top left fig 6: no treatment')
t_final = 120;
time_mesh = linspace(0,t_final,100);
L_t = []; L_d = [];     
I_t = []; I_d = [];      
M_t = []; M_d = [];
L = [time_vec; doses];I = [time_vec; doses];M=[time_vec; doses];
par = mouse();
y0 = [10^6; 5*10^4; 100; 1.1*10^7; 0; 0];
[t,y] = ode45(@(t,y) forwardfunc(t, y, par, L_t, L_d,...
                M_t, M_d, I_t, I_d, time_mesh), time_mesh,y0);
disp('done')
plott(t,y,t_final);
%%
clc
disp('section 5.1: Immune system response to tumor: mouse data')
disp('top right fig 6 : chemotherapy every 14 day to day 90 dose = 1')
t_final = 120;
time_mesh = linspace(0,t_final,100);
L_t = [];         L_d = [];     
I_t = [];         I_d = [];      
M_t=0:14:90;      M_d = ones(1,length(M_t));
par = mouse();
y0 = [10^6; 5*10^4; 100; 1.1*10^7; 0; 0];
[t,y] = ode45(@(t,y) forwardfunc(t, y, par, L_t,L_d,...
                M_t,M_d,I_t,I_d, time_mesh),time_mesh,y0);
disp('done')
plott(t,y,t_final);
%%
clc
disp('section 5.1: Immune system response to tumor: mouse data')
disp('bottom left fig 6 : TIL dose 8*10^8 day 7')
t_final = 120;
time_mesh = linspace(0,t_final,1000);
L_t = [7];      L_d = [8*10^8];     
I_t = [];       I_d = [];      
M_t = [];       M_d = [];
par = mouse();
y0 = [10^6; 5*10^4; 100; 1.1*10^7; 0; 0];
[t,y] = ode45(@(t,y) forwardfunc(t, y, par, L_t, L_d,...
                M_t, M_d, I_t, I_d, time_mesh),time_mesh,y0);
disp('done')
plott(t,y,t_final);
%%
clc
disp('section 5.1: Immune system response to tumor: mouse data')
disp('bottom right fig 6 : combination of TIL 8*10^8 day 7 and chemotherapy every 14 day')
t_final = 120;
time_mesh = linspace(0,t_final,1000);
L_t = [7];      L_d = [8*10^8];     
I_t = [];       I_d = [];      
M_t=0:14:90;    M_d = ones(1,length(M_t));
par = mouse();
y0 = [10^6; 5*10^4; 100; 1.1*10^7; 0; 0];
[t,y] = ode45(@(t,y) forwardfunc(t, y, par, L_t, L_d,...
                M_t, M_d, I_t, I_d, time_mesh),time_mesh,y0);
disp('done')
plott(t,y,t_final);

%%
function f = plott(t,y,t_final)
    y = real(y); ycell = ceil(y(:,1:4)); ycell(ycell == 0) = 1; y(:,1:4)=ycell; 
    subplot(1,2,1)
        semilogy(t,y(:,1:4),'linewidth',2)
        grid on
        legend('Tumor(T)','NK cells(N)','CD8+T cell(L)','circulating lymphocytes(C)')
        xlim([0 t_final])
    subplot(1,2,2)
        plot(t,y(:,5:6),'linewidth',2)
        grid on
        legend('chemotherapy(M)','immunotherapy(I)')
        xlim([0 t_final])
    f = 0;
end