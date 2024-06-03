%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Case B
%%%% Sys-C 
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% simulation  setting 
N      = 400000;
N_eva  = 10000;
T_max  = 100;   

% divergence run
T_index = zeros(1,T_max); 
T_start = 1;
T_end   = 40;
T_num   = (T_end - T_start) + 1;

% Broadband reference noise ==> AR(1) model
c_AR   = 0.75; 
load x_org_ALL_COR.mat; 

% Additive noise in primary noise:v_p(n)
load v_p_ALL_COR.mat;    %std_p    = 0.1;

% AWGN
load v_o_ALL_COR.mat;    % std_vo = 0.50;

% Narrowband noise: p_2(n)
omega   = [0.10   0.15    0.20]*pi;
[wk, q] = size(omega);
Amp     = [1.0   0.5   0.25];
Phs     = [0    pi/6    pi/3];

p_2_org = zeros(1, N);
for i=1:q
    p_2_org = p_2_org + Amp(i) * sin ( omega(i) * (0:N-1) + Phs(i) );
end
    
% Primary path: P(z)
M_p     = 61;
cutoff  = 0.6*pi;
s_p     = fir1(M_p-1, cutoff/pi);

% Secondary path (SP): S(z)
M       = 31;
s_1st   = fir1(M-1, cutoff/pi);
M2      = 33;
s_2nd   = fir1(M2-1, cutoff/pi);

s = zeros( max(M, M2), N);  

for n=1:N
    if n < N/2+1
        s(1:M, n) = s_1st';
    else
        s(1:M2, n) = s_2nd';
    end
end

% W1(z)
L_1  = 71;
mu_1 = 0.000375;

% W2(z)
L_2   = 30;
mu_2  = 0.000375;

% H(z)
L_h   = 71;
mu_h  = 0.000500;

% online secondary path modeling (OSPM)
hat_M  = max(M, M2) + 10;  
mu_s   = 0.000325;     

%Q(z): linear prediction filter 
D     = hat_M + 1;   
L_LP  = q * 40;      
mu_LP = 0.000400;    

% G_s(n): AWGN scaling 
betaa   = 0.9995;    
Gamma   = 1;         
G_bound = 5;         

% refreshment decision, parameters
lambda  = 0.99;                 
T_p     = 50;                   
alphh   = 1.05;                 

% data saving set
save_e_Sys_C      = zeros(1, N);
save_e2_Sys_C     = zeros(1, N);
save_SPMSE_Sys_C  = zeros(1, N);  
save_G_Sys_C      = zeros(1, N); 

save_NRP_Sys_C    = zeros(2, T_num);

for t=T_start:T_end

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
    x_1_org   = zeros(1, N);
    x_1       = zeros(1, N);
    p_1       = zeros(1, N);
    v_p       = v_p_ALL_COR(t, 1:N);
      
    y         = zeros(1, N);
    y_p       = zeros(1, N);
    e         = zeros(1, N);
    
    % OSPM
    v_o   = v_o_ALL_COR(t, 1:N);
    
    y_s    = zeros(1, N);
    e_s    = zeros(1, N);
    hat_s  =  zeros(hat_M, N);
    
    % Q(z): LPF
    y_LP   = zeros(1, N);
    e_LP   = zeros(1, N);
    c_LP   = zeros(L_LP, N);

    %scaler
    v      = zeros(1, N);
    G_B    = zeros(1, N);
    G_N    = zeros(1, N);
    G_BN   = zeros(1, N);
    G      = zeros(1, N);
    
    % Refreshment
    P_e     = zeros(1, N);
    P_e_T   = zeros(1, N/T_p);  
    
    % reference signal x_1(n) with AR(1) model
    x_1_org(1:N) = x_org_ALL_COR(t, 1:N);
    for n=2:N
        x_1(n) = c_AR(1) * x_1(n-1) + x_1_org(n);
    end

    % broadband noise p_1(n)
    for n1=1:N
        for j=0:M_p-1
            if n1-j > 0
                p_1(n1) = p_1(n1) + s_p(j+1) * x_1(n1-j);
            end
        end
    end
    
     % adjustment
    P_rate = sqrt( var(p_1) ) / sqrt( Amp * Amp' / 2 );   
    p_2(1:N) = p_2_org * P_rate;                          
          
      
     %%%%%%%%%%%%%%%%%%%%%%
     % main loop: n loop
     %%%%%%%%%%%%%%%%%%%%%%
     for n=1:N
         % y_1(n)
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
             x_2(n) = x_2(n) + e(n-1)-y_h(n-1)-y_s(n-1);
         end
         
         % y_2(n)
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
         
         % G_s(n)
         if n-1 > 0
             G_BN(n) = betaa * G_BN(n-1) + (1-betaa) * ...
                    ( (abs( y_h(n-1) ))^(Gamma) +  (abs( y_LP(n-1) ))^(Gamma) );
         end
         
         G(n) = G_BN(n);
         
         if G(n) > G_bound
             G(n) = G_bound;
         end


         % v(n)
         v(n) = v_o(n) * G(n);
         
         % y(n)
         y(n) = y_1(n) + y_2(n) - v(n);

         % y_p(n)
         for m=0:max(M, M2)-1
             if n-m > 0
                 y_p(n) = y_p(n) + s(m+1,n) * y(n-m); 
             end
         end
                
         % e(n)
         e(n) = p_1(n) + p_2(n) + v_p(n) - y_p(n);
         
         % e_h(n)
         e_h(n) = e(n) - y_h(n);
         
         % y_LP(n)
         for j=0:L_LP-1
             if n-D-j > 0
                 y_LP(n) = y_LP(n) + c_LP(j+1, n) * e_h(n-D-j);
             end
         end
         
         % e_LP(n)
         e_LP(n) = e_h(n) - y_LP(n);        
         
         % y_s(n)
         for m=0:hat_M-1
             if n-m > 0
                 y_s(n) = y_s(n) + hat_s(m+1, n) * v(n-m);
             end
         end
         
         % error signal of OSPM
         e_s(n) = e_LP(n) - y_s(n);
         
         % Update 
         for m=0:hat_M-1
             if n-m > 0
                 hat_x_1(n) = hat_x_1(n) + hat_s(m+1, n) * x_1(n-m);
                 hat_x_2(n) = hat_x_2(n) + hat_s(m+1, n) * x_2(n-m);
             end
         end
         
         % Update controller W1(z)
         for j=0:L_1-1
             if n-j > 0
                 w1(j+1, n+1) = w1(j+1, n)  + mu_1 * y_h(n) * hat_x_1(n-j);
             end
         end
         
         % Update controller W2(z)
         for j=0:L_2-1
             if n-j > 0
                 w2(j+1, n+1) = w2(j+1, n)  + mu_2 * y_LP(n) * hat_x_2(n-j);
             end
         end
         
         % Update SF H(z)
         for j=0:L_h-1
             if n-j > 0
                 h(j+1, n+1) = h(j+1, n)  + mu_h * e_s(n) * x_1(n-j);
             end
         end
         
         % Update LPF Q(z)
         for j=0:L_LP-1
             if n-D-j > 0
                 c_LP(j+1, n+1) = c_LP(j+1, n)  + mu_LP * e_s(n) * e_h(n-D-j);
             end
         end
          
         % Update OSPM hat{S}(z)
         for m=0:hat_M-1
             if n-m > 0
                 hat_s(m+1,n+1) = hat_s(m+1, n) + mu_s * e_s(n) * v(n-m);
             end
         end
         
         %  refreshment 
         if n-1 >0
             P_e(n) = lambda * P_e(n-1) + (1-lambda) * e(n)^2;
         end
         
         amari = mod(n, T_p);
         if amari == 0 && n/T_p > 1
             P_e_T(n/T_p)  = lambda *  P_e_T( (n/T_p)-1 ) ...
                             + (1-lambda) * sum( P_e(n-T_p+1 : n) );
             % refreshment
             if  n > N/5  && P_e_T(n/T_p) > alphh  * P_e_T( (n/T_p) - 1 )
                 refresh_poi = n;
                 
                   w1(1:L_1, n+1)       = zeros(L_1, 1);  % W1
                   w2(1:L_2, n+1)       = zeros(L_2, 1);  % W2
                     h(1:L_h, n+1)       = zeros(L_h, 1);  % SF (H)
               c_LP(1:L_LP, n+1)    = zeros(L_LP, 1); % LPF;
                 
                 hat_s(1:hat_M, n+1) = zeros(hat_M, 1);   % hat{S}(z)
                 
                 G(n)          = zeros(1, 1);
                 
                 disp( 'Refreshment is done (SPM) ! ');
                 
             end
         end % refreshment end
                

     end   % n loop end

     % index, 
     var_temp = var(e);
     if isnan(var_temp)
         T_index(t) = 1;
     end
 
     % saving
     save_e_Sys_C   = save_e_Sys_C + e(1:N)/T_num;
     save_e2_Sys_C  = save_e2_Sys_C  + e(1:N) .* e(1:N) / T_num;

     for n=1:N
         if hat_M > M
             wk_err_s(1:M)       = hat_s(1:M, n)' - s(1:M,n)';
             wk_err_s(M+1:hat_M) = hat_s(M+1:hat_M, n)' - 0;
         else
             wk_err_s(1:hat_M)   = hat_s(1:hat_M, n)' - s(1:hat_M,n)';
             wk_err_s(hat_M+1:M) = 0 - s(hat_M+1:M,n)';
         end
         
         save_SPMSE_Sys_C(n)    = save_SPMSE_Sys_C(n) + ( wk_err_s * wk_err_s') / T_num;                
     end 

     save_G_Sys_C  = save_G_Sys_C + G(1:N)/T_num;

     save_NRP_Sys_C(1,t) = 10*log10( var(e(N/2-N_eva:N/2)) ...
         / ( var( p_1(N/2-N_eva:N/2) ) + var(p_2(N/2-N_eva:N/2)) + var(v_p(N/2-N_eva:N/2)) ) );

     save_NRP_Sys_C(2,t) = 10*log10( var(e(N-N_eva:N)) ...
         / ( var( p_1(N-N_eva:N) ) + var(p_2(N-N_eva:N)) + var(v_p(N-N_eva:N)) ) );
    
end % t loop end

%saving MSE of  OSPM
for n=1:N
    if n<N/2+1
        save_SPMSE_Sys_C(n)    = save_SPMSE_Sys_C(n)/(s_1st*s_1st');
    else
        save_SPMSE_Sys_C(n)    = save_SPMSE_Sys_C(n)/(s_2nd*s_2nd');
    end   
end

% Ploting
load Case_B_Sys_A.mat;
load Case_B_Sys_B.mat;

% Ploting
N_end  = 400000; 

figure (2)
subplot(4,1,1);
plot(1:N_end, p_1(1:N_end)+p_2(1:N_end)+v_p(1:N_end), '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
%xlabel('p_1(n)+p_2(n)+v_p(n)');

subplot(4,1,2);
plot(1:N_end, save_e_Sys_A(1:N_end) , '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [ 1  N  -4  4]  );
%xlabel('e(n)');

subplot(4,1,3);
plot(1:N_end, save_e_Sys_B(1:N_end) , '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [ 1  N  -4  4]  );
%xlabel('y_h(n)');

subplot(4,1,4);
plot(1:N_end, save_e_Sys_C(1:N_end) , '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
axis( [ 1  N  -4  4]  );
%xlabel('y_{LP}(n)');
xlabel('Iteration number n', ...
       'FontName', 'Arial', 'FontWeight', 'Bold', 'FontSize', 16);

figure (3)
plot(1:N, 10*log10(save_e2_Sys_A(1:N)), '-', ...
     1:N, 10*log10(save_e2_Sys_B(1:N)), '-', ...-
     1:N, 10*log10(save_e2_Sys_C(1:N)), '-');
axis([0 N -25 13 ]);
legend('Sys-A','Sys-B','Sys-C(Proposed)','Location','Northwest');
legend('boxoff')
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
ylabel('Residual noise power [dB]',...
       'FontName', 'Arial', 'FontWeight', 'Bold', 'FontSize', 16);
xlabel('Iteration number n', ...
       'FontName', 'Arial', 'FontWeight', 'Bold', 'FontSize', 16);

figure (4);
subplot(2,1,1)
plot(1:T_num, save_NRP_Sys_C(1,1:T_num), 'o-');
ylabel('NRP(1st half) [dB]')

subplot(2,1,2)
plot(1:T_num, save_NRP_Sys_C(2,1:T_num), 'o-');
ylabel('NRP(2nd half) [dB]')
    
figure (5);
plot(1:N_end, 10*log10( save_SPMSE_Sys_B(1:N_end) ), '-', ...
    1:N_end, 10*log10( save_SPMSE_Sys_C(1:N_end) ), '--');
axis([1  N_end -20  2]);
legend('Sys-B','Sys-C');
legend('boxoff')
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth',1);
xlabel('Iteration number n', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('OSPM MSE [dB]', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);

    
figure (6);
plot(1:N_end, save_G_Sys_C(1:N_end), '-');
set(gca, 'FontSize', 16, 'FontName', 'Arial','LineWidth', 1);
xlabel('Iteration number n', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);
ylabel('Scaling factor G_s(n)', ...
         'fontName', 'Arial', 'fontWeight', 'Bold', 'fontSize', 16);


disp('Mean NRP(1) ==>'); 
disp(mean(save_NRP_Sys_C(1,:)));

disp('Mean NRP(2) ==>'); 
disp(mean(save_NRP_Sys_C(2,:)));




