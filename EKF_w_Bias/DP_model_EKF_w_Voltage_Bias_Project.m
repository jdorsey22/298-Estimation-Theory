
%% Extended Kalman Filter: DP Model with Voltage Bias

clear all, clc
% Import Battery Parameters 
BatteryParams_w_Bias    

Voc0 = 3.435;
% Import Kalman Parameters
KalmanParams_w_Bias

% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')
% load('IV_data_nonlinear.mat')
load('Sim_Truth_FifthOrder_w_Bias.mat')

Q= 1; 
R = 1000;

%% Observability & Controlability Analysis: 

O = obsv(Ad,Cd);

rank(O)

% Observability Matrix is of Full Rank: System is fully Observable  

Cm = ctrb(Ad,Bd);

rank(Cm)

% Controlability Matrix is of Full Rank: System is fully Controlable

%%
% Initial Conditions: 
P1(1) = 0;           % Covariance 
P4(1) = 0;           % Covariance 

x1(1) = 1;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1
x3(1) = 0;          % Vc2
x4(1) = 0 ;          % Vbias

x1_hat(1) = 1; 
x4_hat(1) = 0; 

for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1); % soc
    x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1); % Vc1
    x3(k) = Ad(3,3)*x3(k-1) + Bd(3,1)*I(k-1); % Vc2
    x4(k) = Ad(4,4)*x4(k-1) + Bd(4,1)*I(k-1); % Vbias
    
    
    % Model Prediction: 
    x1_hat_prev = Ad(1,1)*x1_hat(k-1) + Bd(1,1)*I(k-1);
    x4_hat_prev = Ad(4,4)*x4_hat(k-1) + Bd(4,1)*I(k-1); 
    
    C_ek1 = interp1(soc_intpts_OCV_slope', OCV_slope_intpts, x1_hat_prev);
    C_ek4 = 1; 
    

    P_prev1 = A_ek1*P1(k-1)*A_ek1'+ E_ek1*Q*E_ek1';
    P_prev4 = A_ek4*P4(k-1)*A_ek4'+ E_ek4*Z*E_ek4';

    
   % Measurement Update: 
   V_hat(k) = interp1(soc_intpts_OCV',OCV_intpts,x1_hat_prev) - I(k-1)*R0 - x2(k-1) - x3(k-1) + x4_hat_prev;
   % V_hat(k) = Voc0 + alpha*x1_hat_prev - I(k-1)*R0 - x2(k-1)- x3(k-1);
    
    L1 = P_prev1*C_ek1'*inv(C_ek1*P_prev1*C_ek1'+ F_ek1*R*F_ek1');
    L4 = P_prev4*C_ek4'*inv(C_ek4*P_prev4*C_ek4'+ F_ek4*R*F_ek4');

    x1_hat(k) = x1_hat_prev + L1*(V(k)-V_hat(k));
    x4_hat(k) = x4_hat_prev + L4*(V(k)-V_hat(k));

    P1(k) = P_prev1 - L1*C_ek1*P_prev1;
    P4(k) = P_prev4 - L4*C_ek4*P_prev4;
end 



figure();
hold on 
plot(t,SOC_act)
plot(t,x1_hat)
plot(t,x1)
title('Problem #2 Extended Kalman Filter: SOC Results (Jonathan Dorsey)'); 
xlabel('Time (seconds)'); 
ylabel('State of Charge (SOC)'); 

legend('SOC Act','SOC Est','SOC_ OL');


figure(); 
plot(t,x4_hat)

%%

Error = SOC_act-x1_hat';
sigma = sqrt(P(end)); 

int = linspace(-.02,.02, 3500); 

[f,x] = hist(Error,int); 
fx_pdf = pdf('norm',x,0,sigma); 

figure(); hold on; 
plot(x,f);
plot(x,fx_pdf,'r'); 
xlabel('Error') 
ylabel('Prob. Density');
legend('SOC Error','~N(0,P)'), title('SOC Error Distribution (Jonathan Dorsey)');



%% RMS Error Optimation: 

R_list = [.1,1,10,100,1000,2000,5000,10000]; 


for l = 1:length(R_list)
    
    Q= 1; 
R = R_list(l);

% Initial Conditions: 
P(1) = 0;           % Covariance 
x1(1) = 1;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1
x3(1) = 0;          % Vc2

x1_hat(1) = 1; 

for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1); % soc
    x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1); % Vc1
    x3(k) = Ad(3,3)*x3(k-1) + Bd(3,1)*I(k-1); % Vc2
    
    
    % Model Prediction: 
    x1_hat_prev = Ad(1,1)*x1_hat(k-1) + Bd(1,1)*I(k-1);
    
    C_ek = interp1(soc_intpts_OCV_slope', OCV_slope_intpts, x1_hat_prev);

    P_prev = A_ek*P(k-1)*A_ek'+ E_ek*Q*E_ek';
    
   % Measurement Update: 
   V_hat(k) = interp1(soc_intpts_OCV',OCV_intpts,x1_hat_prev) - I(k-1)*R0 - x2(k-1) - x3(k-1);
   % V_hat(k) = Voc0 + alpha*x1_hat_prev - I(k-1)*R0 - x2(k-1)- x3(k-1);
    
L = P_prev*C_ek'*inv(C_ek*P_prev*C_ek'+ F_ek*R*F_ek');
    
    x1_hat(k) = x1_hat_prev + L*(V(k)-V_hat(k));
    P(k) = P_prev - L*C_ek*P_prev;
    
end 

    
error_rms(l) = sqrt(mean((SOC_act-x1_hat).^2)); 
    
    
    
end 




