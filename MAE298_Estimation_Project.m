%% MAE 298 Estimation Project Matlab Script


%% Extended Kalman Filter: DP Model: 

clear all, clc
% Import Battery Parameters 
BatteryParams    

Voc0 = 3.435;
% Import Kalman Parameters
KalmanParams

% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')
load('IV_data_nonlinear.mat')
% load('SimTruth1.mat')

Q= 4; 
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



figure();
hold on 
plot(t,SOC_act)
plot(t,x1_hat)
plot(t,x1)
title('Problem #2 Extended Kalman Filter: SOC Results (Jonathan Dorsey)'); 
xlabel('Time (seconds)'); 
ylabel('State of Charge (SOC)'); 

legend('SOC Act','SOC Est','SOC_ OL');

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

%% Simulated Truth : 
clear all, clc
% Import Battery Parameters 
BatteryParams    

Voc0 = 3.435;
% Import Kalman Parameters
KalmanParams

% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')
load('IV_data_nonlinear.mat')

P(1) = 0;           % Covariance 
x1(1) = 1;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1
x3(1) = 0;          % Vc2

x1_hat(1) = x1(1); 

var = 3.5*10^-7;
% var = 2.5*10^-6;

% var = 0;

for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1)+normrnd(0,sqrt(var)); % soc
    x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1) +normrnd(0,sqrt(var)); % Vc1
    x3(k) = Ad(3,3)*x3(k-1) + Bd(3,1)*I(k-1)+ normrnd(0,sqrt(var)); % Vc2
    
     V_truth(k) = interp1(soc_intpts_OCV',OCV_intpts,x1(k-1)) - I(k-1)*R0 - x2(k-1)- x3(k-1)+normrnd(0,sqrt(R));

end 



figure()
hold on 
plot(t,x1)
plot(t,SOC_act)
legend('Simulated Truth','Lin_SOC_act');

figure(); 
plot(t,V_truth)




