%% Fifth Order Simulation Truth Generation Script: 


%% Set Parameters: 

C1 = 5630; 
C2 = 54277; 
R1 = .00064; 
R2 = .00824;  
R0 = .02402; 
alpha = .65; 
Cbat = 5*3600;  

Tau1 = C1*R1; 
Tau2 = C2*R2; 

dt = .1; 

%% System Dynamics

% Linear State Dynamics: Dual Polarity Model 

% Continuous Time Model: 
A_c = [0       0           0     ; ...
       0  (-1/(R1*C1))     0     ; ... 
       0       0    (-1/(R2*C2))]; 
B_c = [(-1/Cbat); (1/C1); (1/C2)]; 
C_c = [alpha -1 -1];
D_c = [-R0]; 

Ad = [1       0       0      ;
      0  (1-dt/Tau1)  0      ;
      0       0  (1-dt/Tau2)];

 Bd = [-dt/Cbat
       dt/C1
       dt/C2];

% Discrete Time Model: 
% 
% Sys_cont = ss(A_c,B_c,C_c,D_c); 
% 
% Sys_dis = c2d(Sys_cont,dt); 

KalmanParams
% Load Battery Measurements 
load('OCV_table.mat')
load('OCV_slope_table.mat')
load('IV_data_nonlinear.mat')

%% State/Output Simulation with Process/Measurement Noise (Truth) 

P(1) = 0;           % Covariance 
x1(1) = 1;          % SOC - Battery Fully Charged 
x2(1) = 0;          % Vc1
x3(1) = 0;          % Vc2

x1_hat(1) = x1(1); 

var1 = 3.5*10^-7;
var2 = 1*10^-7;
% var4 = 2*10^-6;
var4 = 1.5*10^-7;

var3 = 0;

for k = 2:1:length(t)
    
    x1(k) = Ad(1,1)*x1(k-1) + Bd(1,1)*I(k-1)+ normrnd(0,sqrt(var2)); % soc
    x2(k) = Ad(2,2)*x2(k-1) + Bd(2,1)*I(k-1) +normrnd(0,sqrt(var4)); % Vc1
    x3(k) = Ad(3,3)*x3(k-1) + Bd(3,1)*I(k-1)+ normrnd(0,sqrt(var4)); % Vc2
    
    V_truth(k) = interp1(soc_intpts_OCV',OCV_intpts,x1(k-1))-I(k-1)*R0-x2(k-1)+normrnd(0,sqrt(R));
end 

figure()
hold on 
plot(t,x1)
plot(t,SOC_act)
legend('Simulated Truth','Lin_SOC_act');

figure(); 
plot(t,V_truth)

%% NaN Determination for Preventing Interpolation Fuckups

thing = isnan(V_truth);

counter =0; 
for k=1:length(t)
    if thing(k)==1
%        display('yes') 
       
       counter = counter +1;
       
    end     
end 
counter
SOC_act(thing) = 0;

%% Relabel New truth Data

clear V SOC_act 

V = V_truth; 
SOC_act = x1; 

save('ThirdOrder_Truth_BestNAN','V','SOC_act','t','I'); 


%% 
hold on 
plot(t,SOC_act)