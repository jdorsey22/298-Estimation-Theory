%% Extended Kalman Filter Parameters:

wk_mean = 0; 
Q = 2.5*10^-7;
vk_mean = 0; 
R = 1*10^-4;
Z = R; 


%% X1 Hat Kalman Parameters 
A_ek1 = 1 ;
E_ek1 = 1; 
F_ek1 = 1; 

%% X4 Hat Kalman Parameters 

A_ek4 = 1; 
E_ek4 = 1; 
F_ek4 = 1; 

%%



Ak = 1;
Bk = -dt/Cbat;

