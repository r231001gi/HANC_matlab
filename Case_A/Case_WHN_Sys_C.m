%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HANC Sys-C  - Case A
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N     = 100000;
N_eva = 10000;
T     = 40;

% Broadband primary noise
load x_org_ALL_WHN.mat; 

% Additive noise in primary noise:v_p(n)
load v_p_ALL.mat;    %std_p       = 0.1;

% AWGN
load v_o_ALL.mat;    % std_vo = 0.50;

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
M       = 31;
s       = fir1(M-1, cutoff/pi);

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
mu_s   = 0.00725;  

%Q(z)  - LPF
D     = hat_M + 1;
L_LP  = q * 40;
mu_LP = 0.000425;

betaa  = 0.9995;   % [0, 1)
Gamma  = 1;        %  {1, 2, 3, 4}

save_e_Sys_C     = zeros(1, N);
save_e2_Sys_C    = zeros(1, N);
save_SPMSE_Sys_C = zeros(1, N);
save_G_Sys_C     = zeros(1, N);

save_NRP_Sys_C  = zeros(1, T);

for t=1:T
    disp('t ===>'); disp(t);
    
    % W1
    y_1       = zeros(1, N);
    hat_x_1   = zeros(1, N);
    w1        = zeros(L_1, N);
    
    % W2
    y_2       = zeros(1, N);
    x_2       = zeros(1, N);
    hat_x_2   = zeros(1, N);
    w2        = zeros(L_2, N);
    
    % H(z)
    y_h       = zeros(1, N);
    e_h       = zeros(1, N);
    h         = zeros(L_h, N);
    
    % HANC
    p_1       = zeros(1, N);
    v_p       = v_p_ALL(t, 1:N);
    
    x_1       = zeros(1, N);
     
    y         = zeros(1, N);
    y_p       = zeros(1, N);
    e         = zeros(1, N);
    
    % OSPM
    v_o    = v_o_ALL(t, 1:N);
    
    y_s    = zeros(1, N);
    e_s    = zeros(1, N);
    hat_s  =  zeros(hat_M, N);
    
    % Q(z) - LPF
    y_LP   = zeros(1, N);
    e_LP   = zeros(1, N);
    c_LP   = zeros(L_LP, N);

    %gain
    v      = zeros(1, N);
    G_B    = zeros(1, N);
    G_N    = zeros(1, N);
    G      = zeros(1, N);
    G_BN   = zeros(1,N);


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
             x_2(n) = x_2(n) + e(n-1)- y_h(n-1)- y_s(n-1);
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
         
         % AWGN scaling
         if n-1 > 0
             G_BN(n) = betaa * G_BN(n-1) + (1-betaa) * ...
                    ( (abs( y_h(n-1) ))^(Gamma) +  (abs( y_LP(n-1) ))^(Gamma) );
         end

         G(n) = G_BN(n);
         
         % v(n)
         v(n) = v_o(n) * G(n);
         
         % y(n), y_p(n)
         y(n) = y_1(n) + y_2(n) - v(n);
         for m=0:M-1
             if n-m > 0
                 y_p(n) = y_p(n) + s(m+1) * y(n-m);
             end
         end
                
         % e(n)
         e(n) = p_1(n) + p_2(n) + v_p(n) - y_p(n);
         
         % e_h(n)
         e_h(n) = e(n) - y_h(n);
         
         % Q(z) - LPF
         for j=0:L_LP-1
             if n-D-j > 0
                 y_LP(n) = y_LP(n) + c_LP(j+1, n) * e_h(n-D-j);
             end
         end
         
         e_LP(n) = e_h(n) - y_LP(n);        

         % OSPM
         for m=0:hat_M-1
             if n-m > 0
                 y_s(n) = y_s(n) + hat_s(m+1, n) * v(n-m);
             end
         end
         
         e_s(n) = e_LP(n) - y_s(n);

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
                 w2(j+1, n+1) = w2(j+1, n)  + mu_2 * y_LP(n) * hat_x_2(n-j);
             end
         end
         
         % Update H
         for j=0:L_h-1
             if n-j > 0
                 h(j+1, n+1) = h(j+1, n)  + mu_h * e_s(n) * x_1(n-j);
             end
         end
         
         % Update Q(z) - LPF
         for j=0:L_LP-1
             if n-D-j > 0
                 c_LP(j+1, n+1) = c_LP(j+1, n)  + mu_LP * e_s(n) * e_h(n-D-j);
             end
         end         
         
         % Update hat{S}(z)
         for m=0:hat_M-1
             if n-m > 0
                 hat_s(m+1,n+1) = hat_s(m+1, n) + mu_s * e_s(n) * v(n-m);
             end
         end

     end   % n loop end
     
     % 
     save_e_Sys_C    = save_e_Sys_C + e(1:N)/T;
     save_e2_Sys_C   = save_e2_Sys_C  + e(1:N) .* e(1:N) / T;
     
     for n=1:N
         if hat_M > M
             wk_err_s(1:M)             = hat_s(1:M, n)' - s(1:M);
             wk_err_s(M+1:hat_M) = hat_s(M+1:hat_M, n)' - 0;
         else
             wk_err_s(1:hat_M)      = hat_s(1:hat_M, n)' - s(1:hat_M);
             wk_err_s(hat_M+1:M) = 0 - s(hat_M+1:M);
         end
         
         save_SPMSE_Sys_C(n)    = save_SPMSE_Sys_C(n) + ( wk_err_s * wk_err_s') / T;                
     end 
     
     save_G_Sys_C  = save_G_Sys_C + G(1:N)/T;

     save_NRP_Sys_C(t) = 10*log10( var(e(N-N_eva:N)) ...
         / ( var( p_1(N-N_eva:N) ) + var(p_2(N-N_eva:N)) + var(v_p(N-N_eva:N)) ) );
    
end % t loop end

% Ploting
load Case_WHN_Sys_A.mat;
load Case_WHN_Sys_B.mat;

% Ploting
N_end  = 50000;

figure (2);
subplot(4,1,1);
plot(1:N_end, p_1(1:N_end)+p_2(1:N_end)+v_p(1:N_end), '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [0   N_end  -5  5] );

subplot(4,1,2);
plot(1:N_end, save_e_Sys_A(1:N_end) , '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [0   N_end  -2  2] );

subplot(4,1,3);
plot(1:N_end, save_e_Sys_B(1:N_end) , '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [0   N_end  -2  2] );

subplot(4,1,4);
plot(1:N_end, save_e_Sys_C(1:N_end) , '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [0   N_end  -2  2] );
xlabel('Iteration number n', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);


figure (3)
plot(1:N_end, 10*log10(save_e2_Sys_A(1:N_end)), '-', ...
     1:N_end, 10*log10(save_e2_Sys_B(1:N_end)), '-', ...
     1:N_end, 10*log10(save_e2_Sys_C(1:N_end)), '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
ylabel('Residual noise power [dB]', ...
     'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
xlabel('Iteration number n', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

figure (4);
plot(1:T, save_NRP_Sys_C(1:T), 'o-');
ylabel('NRP [dB]')


figure (5);
plot(1:N_end, 10*log10( save_SPMSE_Sys_B(1:N_end)/(s * s') ), '-', ...
     1:N_end, 10*log10( save_SPMSE_Sys_C(1:N_end)/(s * s') ), '--');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('OSPM MSE [dB]', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
    
figure (6)
plot(1:N_end, save_G_Sys_C(1:N_end), '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('Scaling factor G_s(n)', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);     
     
disp('mean NRP');
disp(mean(save_NRP_Sys_C));


