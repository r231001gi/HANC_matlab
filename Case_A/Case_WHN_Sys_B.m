%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HANC Sys-B  - Case A
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N      = 100000;
N_eva  = 10000;
T      = 40;

% Broadband primary noise
load x_org_ALL_WHN.mat; 

% Additive noise in primary noise:v_p(n)
load v_p_ALL.mat;    %std_p       = 0.1;

% AWGN
std_vo  = 0.50;
FLG_awgn = 0;
if FLG_awgn == 0
    v_o_ALL = randn(100, 200000) * std_vo;
    save v_o_ALL.mat   v_o_ALL;
else
    load v_o_ALL.mat; 
end

% Narrowband noise p_2(n)
omega     = [0.10   0.15    0.20]*pi;
[wk, q]   = size(omega);
Amp       = [1.0   0.5   0.25];
Phs       = [0    pi/6    pi/3];

p_2_org = zeros(1, N);
for i=1:q
    p_2_org = p_2_org + Amp(i) * sin ( omega(i) * (0:N-1) + Phs(i) );
end
    
% Primary path P(z)
M_p     = 61;
cutoff  = 0.6*pi;
s_p     = fir1(M_p-1, cutoff/pi);

% Secondary path S(z)
M     = 31;
s     = fir1(M-1, cutoff/pi);

% W1(z)
L_1  = 71;
mu_1 = 0.00075;

%W2(z)
L_2   = 30;
mu_2  = 0.00075;

%H(z)
L_h   = 71;
mu_h  = 0.0010;

% OSPM
hat_M  = M + 10;
mu_s   = 0.0005;

save_e_Sys_B      = zeros(1, N);
save_e2_Sys_B     = zeros(1, N);
save_SPMSE_Sys_B  = zeros(1, N);

save_NRP_Sys_B  = zeros(1, T);

for t=1:T
    disp('t ===>'); disp(t);
    
    % W1
    y_1      = zeros(1, N);
    hat_x_1  = zeros(1, N);
    w1       = zeros(L_1, N);
    
    % W2
    y_2      = zeros(1, N);
    x_2      = zeros(1, N);
    hat_x_2  = zeros(1, N);
    w2       = zeros(L_2, N);
    
    % H(z)
    y_h   = zeros(1, N);
    e_h   = zeros(1, N);
    h     = zeros(L_h, N);
    
    % HANC
    p_1   = zeros(1, N);   
    v_p   = v_p_ALL(t, 1:N);
    
    x_1  = zeros(1, N);   
    y    = zeros(1, N);
    y_p  = zeros(1, N);
    e    = zeros(1, N);
    
    % OSPM
    v_o    = v_o_ALL(t, 1:N);
    
    y_s    = zeros(1, N);
    e_s    = zeros(1, N);
    hat_s  =  zeros(hat_M, N);

    % adjustment
    x_1(1:N) = x_org_ALL(t, 1:N);
    for n1=1:N
        for j=0:M_p-1
            if n1-j > 0
                p_1(n1) = p_1(n1) + s_p(j+1) * x_1(n1-j);
            end
        end
    end

    P_rate =  sqrt( var(p_1) ) / sqrt( Amp * Amp' / 2 );   
    
    p_2(1:N) = p_2_org * P_rate;         
    
           
     % main loop
     for n=1:N
         % W1
         for j=0:L_1-1
             if n-j > 0
                 y_1(n) = y_1(n) + w1(j+1, n) * x_1(n-j);
             end
         end
         
         % x_2(n)
         for m=0:hat_M-1
             if n-1-m>0
                 x_2(n) = x_2(n) + hat_s(m+1, n) * y_2(n-1-m);
             end
         end
         
         if n-1>0
             x_2(n) = x_2(n) + e_h(n-1);
         end
         
         % W2
         for j=0:L_2-1
             if n-j > 0
                 y_2(n) = y_2(n) + w2(j+1, n) * x_2(n-j);
             end
         end
         
         % y_h(n)
         for j=0:L_h-1
             if n-j > 0
                 y_h(n) = y_h(n) + h(j+1, n) * x_1(n-j);
             end
         end
         
         % y(n), y_p(n)
         y(n) = y_1(n) + y_2(n) - v_o(n);
         for m=0:M-1
             if n-m > 0
                 y_p(n) = y_p(n) + s(m+1) * y(n-m);
             end
         end
                
         % e(n)
         e(n) = p_1(n) + p_2(n) + v_p(n) - y_p(n);
         
         % e_h(n)
         e_h(n) = e(n) - y_h(n);
         
         % OSPM
         for m=0:hat_M-1
             if n-m > 0
                 y_s(n) = y_s(n) + hat_s(m+1, n) * v_o(n-m);
             end
         end
         
         e_s(n) = e(n) - y_s(n);
         
         % Update preparations
         for m=0:hat_M-1
             if n-m > 0
                 hat_x_1(n) = hat_x_1(n) + hat_s(m+1, n) * x_1(n-m);
                 hat_x_2(n) = hat_x_2(n) + hat_s(m+1, n) * x_2(n-m);
             end
         end
         
         % Update W1
         for j=0:L_1-1
             if n-j > 0
                 w1(j+1, n+1) = w1(j+1, n)  + mu_1 * y_h(n) * hat_x_1(n-j);
             end
         end
         
         % Update W2
         for j=0:L_2-1
             if n-j > 0
                 w2(j+1, n+1) = w2(j+1, n)  + mu_2 * e_h(n) * hat_x_2(n-j);
             end
         end
         
         % Update H
         for j=0:L_h-1
             if n-j > 0
                 h(j+1, n+1) = h(j+1, n)  + mu_h * e_h(n) * x_1(n-j);
             end
         end

         % Update hat{S}(z)
         for m=0:hat_M-1
             if n-m > 0
                 hat_s(m+1,n+1) = hat_s(m+1, n) + mu_s * e_s(n) * v_o(n-m);
             end
         end

         
     end   % n loop end
     
     %
     save_e_Sys_B    = save_e_Sys_B + e(1:N)/T;
     save_e2_Sys_B   = save_e2_Sys_B  + e(1:N) .* e(1:N) / T;
     
     for n=1:N
         if hat_M > M
             wk_err_s(1:M)             = hat_s(1:M, n)' - s(1:M);
             wk_err_s(M+1:hat_M) = hat_s(M+1:hat_M, n)' - 0;
         else
             wk_err_s(1:hat_M)      = hat_s(1:hat_M, n)' - s(1:hat_M);
             wk_err_s(hat_M+1:M) = 0 - s(hat_M+1:M);
         end

         save_SPMSE_Sys_B(n)    = save_SPMSE_Sys_B(n) + ( wk_err_s * wk_err_s') / T;
     end


     save_NRP_Sys_B(t) = 10*log10( var(e(N-N_eva:N)) ...
         / ( var( p_1(N-N_eva:N) ) + var(p_2(N-N_eva:N)) + var(v_p(N-N_eva:N)) ) );
    
     
     
end % t loop end

save Case_WHN_Sys_B.mat  save_e_Sys_B   save_e2_Sys_B  save_SPMSE_Sys_B;

% Ploting
load Case_WHN_Sys_A.mat;

% Ploting
figure (2)
subplot(2,1,1);
plot(1:N, p_1(1:N)+p_2(1:N)+v_p(1:N), '-');
xlabel('p_1(n)+p_2(n)+v_p(n)');

subplot(2,1,2);
plot(1:N, save_e_Sys_B(1:N) , '-');
xlabel('e(n)');


figure(3);
plot(1:N, 10*log10(save_e2_Sys_A(1:N)), '-', ...
     1:N, 10*log10(save_e2_Sys_B(1:N)), '--');
ylabel('Residual noise power [dB]')
axis( [0   N_end  -25  10] );


figure (4);
plot(1:T, save_NRP_Sys_B(1:T), 'o-');
ylabel('NRP [dB]')

figure (5);
plot(1:N, 10*log10( save_SPMSE_Sys_B(1:N)/(s * s') ), '-');

disp('mean NRP');
disp(mean(save_NRP_Sys_B));


