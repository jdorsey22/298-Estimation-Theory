%% Battery Parameters: 






%%

% R1 = 2.85-10^-2; 
% R2 = 4.44*10^-2;
% C1 = 4.78*10^2;
% C2 = 1.83*10^4;
% R0 = 5.55*10^-2;
% 








%%

% C1 = 5630; 
% C2 = 54277; 
% R1 = .00064; 
% R2 = .00824; 
%  R0 = .02402; 
% alpha = .65; 
% Cbat = 5*3600;  

C1 = 2400; 
C2 = 2400; 
R1 = .015; 
R2 = .0015; 
 R0 = .01; 
alpha = .65; 
Cbat = 5*3600; 

Tau1 = C1*R1; 
Tau2 = C2*R2; 

dt = .1; 

%% System Dynamics

% Linear State Dynamics: Dual Polarity Model 

% Continuous Time Model: 
A_c = [0       0         0 ; ...
     0  (-1/(R1*C1))   0 ;... 
     0       0    (-1/(R2*C2))]; 
B_c = [(-1/Cbat); (1/C1); (1/C2)]; 
C_c = [alpha -1 -1 ];
D_c = [-R0]; 

% Discrete Time Model: 

Ad = [1      0        0 ; ...
     0 exp(-dt/Tau1) 0 ; ...
     0      0   exp(-dt/Tau2)]; 
Bd = [(-dt/Cbat); (R1)*(1-exp(-dt/Tau1)); (R2)*(1-exp(-dt/Tau2))]; 
Cd = C_c; 
Dd = D_c; 








