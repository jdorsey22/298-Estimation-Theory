%% Extended Kalman Filter Parameters:

wk_mean = 0; 
Q = 2.5*10^-7;
vk_mean = 0; 
R = 1*10^-4;

A_ek = 1 ;
% C_ek: In main script
E_ek = 1; 
F_ek = 1; 

Ak = 1;
Bk = -dt/Cbat;

