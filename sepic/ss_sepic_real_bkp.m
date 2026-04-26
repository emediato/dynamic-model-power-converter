
%%
%bode(Gntf)
%%
%% REAL VALUES 
Vi_val = 50  ; % Input voltage [V]
Vi=Vi_val;
Vo_val = 100 ; % Output voltage [Vo]
Vo=Vo_val;

D_val = Vo_val / (Vi_val + Vo_val); % SEPIC duty cycle formula

D=D_val;d=D;
fs = 50e3; % [kHz]  1/Ts 
Ts = (1/fs) ; % secs
Vs  = Vi_val*D_val/(1-D_val);
Po = 250 ; % W

Io = Po/Vo ;
Iin = Po/Vi_val ;
IL1_val = Iin; I_l1 = IL1_val;

IL2_val = Io;I_l2= IL2_val;
Rload = Vo/Io;       % load resistance (Ohm)
Ro = Rload ;

delta_Vo = 5/100 ; % 5% 
delta_Q = C * delta_Vo;
Vrms = sqrt(Vi);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

delta_I_L_p = 0.2 ; % 20 % percentage

delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;

L2 = (Vi*D)/(delta_I_l2*fs) ;
L1 = (Vo*(1-D))/(delta_I_l1*fs) ;

delta_V_C1_p = 20/100 ; % 20 percentage
delta_V_C1 = Vi * delta_V_C1_p; 
delta_V_C2 = Vo * delta_Vo; 

C1_val = (Io * D) / (delta_V_C1 * fs);       % Capacitor C1 (F)
C2_val = (Vo * (1-D)) / (8 * delta_V_C2 * fs^2 * L2)       % Capacitor C2 (F)
C1=C1_val;C2=C2_val;

delta_Q_bkp = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;

%% dimensions

%Xn = subs(X, [d, Ro, Vi, rc, rL, rs], [Dn, Rn, Vin, rcn, rLn, rsn])
%VL1avg = Vin;
%VL2avg = Xn(4) - Xn(3);
%IC1avg = abs(Xn(2));

%L1n = VL1avg * Dn / (Xn(1) * delta_iL_p * fs);
%L2n = VL2avg * Dn / (Xn(2) * delta_iL_p * fs);
%C1n = IC1avg * Dn / (Xn(3) * delta_iL_p * fs);
%C2n = abs(Xn(2)) * delta_I_l2 / (4 * Xn(4) * delta_V_C2 * fs);
