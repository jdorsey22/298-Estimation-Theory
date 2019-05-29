%% Battery Model - Felipe Valdez
clear all
close all

%Pre-locating arrays
Vc = zeros(1,length(t));    
SOC = zeros(1,length(t));      
P = zeros(1,length(t));
SOC_ol = zeros(1,length(t)); 

%Battery parameters
C_bat =5*3600; 
R0 = 0.01; %Ohms
Rc = 0.015; %Ohms
Cc = 2400; %F
alpha = 0.65; %V
Vocv0 = 3.435; %V  

%%%%%%%%%%%%%%% Current Input Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load IV_data_nonlinear
load OCV_table
load OCV_slope_table

%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vc(1) = 0;
SOC(1) = 1;
OCV2(1) = Vocv0 + alpha*SOC(1); 
V(1) = OCV2(1); %at t=0, V(t)=OVC(t) since there is no V drop
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
Ak = 1;  Bk = -dt/C_bat; Dk = -R0; Ek = 1; Fk = 1;
Vc(k) = (1-(dt/(Rc*Cc)))*Vc(k-1) + (dt/Cc)*I(k-1);

%Open Loop
SOC_ol(k) = Ak*SOC_ol(k-1)+Bk*I(k);

%%%%%%%%%%%%%%%%%%%%% MODEL PREDICTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
SOC_prev = Ak*SOC(k-1) + Bk*I(k-1);  %Step 1a: State estimate time update 
P_prev = Ak*P(k-1)*Ak'+ Ek*Q*Ek';    %Step 1b: Error covariance time update 

OCV2(k) = interp1(soc_intpts_OCV, OCV_intpts, SOC_prev); %interpolate OCV
V_est = OCV2(k) - Vc(k) - R0*I(k); %Step 1c: Estimate system output 

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
figure, plot(t,e), xlabel('time(s)'), ylabel('error')