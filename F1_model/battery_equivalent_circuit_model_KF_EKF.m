%% KF Battery Model - Felipe Valdez
clear; clc; close all; tic

%%%%%%%%%%%%%%% Current Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load IV_data_linear

%Battery parameters
C_bat =5*3600; 
Rs = 0.01; %Ohms
R1 = 0.015; %Ohms
R2 = 0.015; %Ohms   %decreasing this voltage gives a better estimation...
C1 = 2400; %F  %decreasing C below 2000 provides shitty estimation since the Capacitor deals with nonliniearities 
C2 = 2400; %F
alpha = 0.65; %V
Vocv0 = 3.435; %V

%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%
Vc1 = zeros(1,length(t));       %prelocate array
SOC = zeros(1,length(t));     
P = zeros(1,length(t));
SOC_ol = zeros(1,length(t));   
fx = zeros(1,length(t));       

Vc1(1) = 0;
Vc2(1) = 0;
SOC(1) = 1;
OCV(1) = Vocv0 + alpha*SOC(1); 
V(1) = OCV(1); %at t=0, V(t)=OVC(t) since there is no V drop
mu = 0;                        %zero mean
Q = 2.5e-7;                    %system noise covariance
R = 10^-4;                     %measurement(voltage) noise covariance
SOC_ol(1) = SOC(1);
dt = 0.1;                      %sampling period

%%%%%%%%%%%%%%%%%%%%%%%%% Kalman Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:length(t)
% w = normrnd(mu,Q);
% v = normrnd(mu,R);

%State space 
A = 1;  B = -dt/C_bat;  C = alpha;  D = -Rs;

%Open Loop
SOC_ol(k) = A*SOC_ol(k-1)+B*I(k);

Vc1(k) = (1-(dt/(R1*C1)))*Vc1(k-1) + (dt/C1)*I(k-1);
Vc2(k) = (1-(dt/(R2*C2)))*Vc2(k-1) + (dt/C2)*I(k-1);

%%%%%%%%%%%%%%%%%%%%% MODEL PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
SOC_prev = A*SOC(k-1) + B*I(k-1); %Step 1a: State estimate time update
P_prev = A*P(k-1)*A'+ Q;  %Step 1b: Error covariance time update 
V_est = alpha*SOC_prev - Vc1(k) - Vc2(k) - Rs*I(k) + Vocv0; %Step 1c: Estimate system output 

%%%%%%%%%%%%%%%%%%%%% MEASUREMENT UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = P_prev*C'*inv(C*P_prev*C'+R); %Step 2a: Compute Kalman gain 
SOC(k) = SOC_prev + L*(V(k) - V_est); %Step 2b: State estimate measurement update 
P(k) = P_prev - L*C*P_prev; %Step 2c: Error covariance measurement update 

end

figure; plot(t,SOC,t,SOC_act,t,SOC_ol), legend('estimated', 'actual', 'open loop')
xlabel('time'), ylabel('SOC'), title('Felipe Valdez')

%---------------------------------Part c-----------------------------------
syms P_inf

%Algebraic Riccatti equation
R_eqn = A*P_inf*A.' + Q - A*P_inf*C.'*(C*P_inf*C.'+R)^-1*C*P_inf*A.'-P_inf;
sol = solve(R_eqn==0,P_inf); %solve Riccatti eqn
P_inf = eval(sol(2)); %evaluate Riccatti eqn and get solution 2

P_inf_prev = A*P_inf*A.' + Q - A*P_inf*C.'*(C*P_inf*C.'+R)^-1*C*P_inf*A.';

%compare P_inf Riccatti with KF 
Pk_inf = P_inf_prev - P_inf_prev*C.'*(C*P_inf_prev*C.'+R)^-1*C*P_inf_prev
Pk_inf_KF = P(end)

%--------------------------------Part d------------------------------------
e = SOC_act - SOC';
sigma = sqrt(P(end));

interval = linspace(-0.02,0.02,3500);
[f,x] = hist(e,interval); %use hist function and get unnormalized values
fx_pdf = pdf('norm',x,0,sigma); 

figure; plot(x,f,x,fx_pdf,'r');
xlabel('error'), ylabel('probability density')
legend('SOC estimation error','N(0,P)'), title('Felipe Valdez')
figure, plot(t,e), xlabel('time(s)'), ylabel('error'), title('Felipe Valdez')

%% EKF Battery Model 
clear; 

%%%%%%%%%%%%%%% Current Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load IV_data_nonlinear
load OCV_table
load OCV_slope_table

%Pre-locating arrays
Vc = zeros(1,length(t));    
SOC = zeros(1,length(t));      
P = zeros(1,length(t));
SOC_ol = zeros(1,length(t)); 

%Battery parameters
C_bat =5*3600; 
Rs = 0.01; %Ohms
R1 = 0.015; %Ohms
R2 = 0.015; %Ohms   %decreasing this voltage gives a better estimation...
C1 = 2000; %F
C2 = 2400; %F
alpha = 0.65; %V
Vocv0 = 3.435; %V  

%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vc(1) = 0;
Vc2(1) = 0;
SOC(1) = 1;
OCV1(1) = Vocv0 + alpha*SOC(1); 
V(1) = OCV1(1); %at t=0, V(t)=OVC(t) since there is no V drop
mu = 0;                        %zero mean
Q = 2.5e-7;                    %system noise covariance
R = 10^-4;                     %measurement(voltage) noise covariance
SOC_ol(1) = SOC(1);
dt = 0.1;  % sampling period

%%%%%%%%%%%%%%%%%%%%% Extended Kalman Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 2:length(t)
% w = normrnd(mu,Q);
% v = normrnd(mu,R);

%Set State space matrices
Ak = 1;  Bk = -dt/C_bat; Dk = -Rs; Ek = 1; Fk = 1;
Vc(k) = (1-(dt/(R1*C1)))*Vc(k-1) + (dt/C1)*I(k-1);
Vc2(k) = (1-(dt/(R2*C2)))*Vc2(k-1) + (dt/C2)*I(k-1);

%Open Loop
SOC_ol(k) = Ak*SOC_ol(k-1)+Bk*I(k);

%%%%%%%%%%%%%%%%%%%%% MODEL PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
SOC_prev = Ak*SOC(k-1) + Bk*I(k-1);  %Step 1a: State estimate time update 
P_prev = Ak*P(k-1)*Ak'+ Ek*Q*Ek';    %Step 1b: Error covariance time update 

OCV1(k) = interp1(soc_intpts_OCV, OCV_intpts, SOC_prev); %interpolate OCV

V_est = OCV1(k) - Vc(k) - Vc2(k)- Rs*I(k); %Step 1c: Estimate system output 

%%%%%%%%%%%%%%%%%%%%% MEASUREMENT UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ck = interp1(soc_intpts_OCV_slope, OCV_slope_intpts, SOC_prev); %interpolate OCV slope
L = P_prev*Ck'*inv(Ck*P_prev*Ck' + Fk*R*Fk'); %Step 2a: Compute Kalman gain 
SOC(k) = SOC_prev + L*(V(k) - V_est); %Step 2b: State estimate measurement update 
P(k) = P_prev - L*Ck*P_prev; %Step 2c: Error covariance measurement update 

end

%Plots
figure; plot(t,SOC,t,SOC_act,t,SOC_ol), legend('estimated', 'actual', 'open loop')
xlabel('time'), ylabel('SOC'), title('Felipe Valdez')

%------------------------------Part c--------------------------------------

e = SOC_act - SOC';   %SOC error
sigma = sqrt(P(end));

interval = linspace(-0.02,0.02,3500);
[f,x] = hist(e,interval); %use hist function and get unnormalized values
fx_pdf = pdf('norm',x,0,sigma); %get normal dist. pdf

figure; plot(x,f,x,fx_pdf,'r');
xlabel('error'), ylabel('probability density')
legend('SOC estimation error','N(0,P)'), title('Felipe Valdez')
figure, plot(t,e), xlabel('time(s)'), ylabel('error'), title('Felipe Valdez')

