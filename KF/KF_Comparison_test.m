%% 3rd Order Truth & 2nd Order EKF Model Validation

clear all, clc
% Import Battery Parameters 
addpath('C:\Users\felip\Documents\298-Estimation-Theory\Model Validation\DataFiles')

C1 = 1000; 
C2 = 2500; 
% C1 = 2400; 
% C2 = 2400; 
R1 = .015; 
R2 = .0015; 
R0 = .02402; 
alpha = .65; 
Cbat = 5*3600; 
Vocv0 = 3.435; %V


Tau1 = C1*R1; 
Tau2 = C2*R2; 

dt = .1; 

% Continuous Time Model: 
A_c = [0       0         0 ; ...
       0  (-1/(R1*C1))   0 ;... 
       0       0    (-1/(R2*C2))]; 
B_c = [(-1/Cbat); (1/C1); (1/C2)]; 
C_c = [alpha -1 -1 ];
D_c = [-R0]; 

sys = ss(A_c,B_c,C_c,D_c); 

sys_d= c2d(sys,dt);

% Discrete Time Model: 

Ad = [1      0        0 ; ...
     0 exp(-dt/Tau1) 0 ; ...
     0      0   exp(-dt/Tau2)]; 
Bd = [(-dt/Cbat); (R1)*(1-exp(-dt/Tau1)); (R2)*(1-exp(-dt/Tau2))]; 
Cd = C_c; 
Dd = D_c; 

% Import Kalman Parameters
% KalmanParams


wk_mean = 0; 
Q = 2.5*10^-7;

vk_mean = 0; 
R = 1*10^-4;

A_ek = 1 ;
E_ek = 1; 
F_ek = 1; 

% 
% Q = 1; 
% R = 10000;

% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')
load('ThreeRCModel_Validation_Data.mat')

% Initial Conditions: 
P(1) = 0;           % Covariance 
x1(1) = .98;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1
x3(1) = 0;          % Vc2

x1_hat(1) = x1(1); 

for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1); % soc
%     x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1); % Vc1
%     x3(k) = Ad(3,3)*x3(k-1) + Bd(3,1)*I(k-1); % Vc2
    x2(k) = (1-(dt/(R1*C1)))*x2(k-1) + (dt/C1)*I(k-1); % Vc1
    x3(k) = (1-(dt/(R2*C2)))*x3(k-1) + (dt/C2)*I(k-1); % Vc2
    
    % Model Prediction: 
    x1_hat_prev = Ad(1,1)*x1_hat(k-1) + Bd(1,1)*I(k-1);
    
    if(x1_hat_prev >1)
        x1_hat_prev = 1; 
    end 
%     
    C_ek = alpha;

    P_prev = Ad(1,1)*P(k-1)*Ad(1,1)'+ Q;
    
   % Measurement Update: 
   V_hat(k) = alpha*x1_hat_prev - I(k)*R0 - x2(k) - x3(k) + Vocv0;
    
   L = P_prev*Cd(1,1)'*inv(Cd(1,1)*P_prev*Cd(1,1)'+ R');
    
    x1_hat(k) = x1_hat_prev + L*(V(k)-V_hat(k));
    P(k) = P_prev - L*Cd(1,1)*P_prev;
    
end 



figure(1);
hold on 
plot(t,SOC_act)
plot(t,x1_hat)
plot(t,x1)
title('Linear Kalman Filter SOC Estimation'); %KEEP this Loaded FILE
xlabel('Time (seconds)'); 
ylabel('State of Charge (SOC)'); 
legend('SOC_{act}','SOC_{est}','SOC_{OL}'), grid on;

figure(2)
plot(t, V, t ,V_hat), legend('V_{act}', 'V_{est}')
xlabel('Time (seconds)'), ylabel('V (Volts)'), grid on;
title('Linear Kalman Filter SOC Estimation'); 

%%
e = SOC_act - x1_hat';
sigma = sqrt(P(end));

interval = linspace(-0.02,0.02,3500);
[f,x]=hist(e,interval); %use hist function and get unnormalized values
fx_pdf = pdf('norm',x,0,sigma); 

figure; plot(x,f,x,fx_pdf,'r');
xlabel('error'), ylabel('probability density')
legend('SOC estimation error','N(0,P)'), title('Felipe Valdez')
figure, plot(t,e), xlabel('time(s)'), ylabel('error')

