Vin = 200  ; % Input voltage [V]
Vopk = 311 ; % Output voltage [Vo]
fs = 50e3 ; % 1/Ts [kHz] 
Ts  = 1 / fs;
Po = 250 ; % W

delta_I_L_p = 0.2 ; % 20% percentage
delta_V_C_p = .20 ; % 20 % percentage

Dn = 1/2 * (sqrt(4*Vin^2+Vopk^2) - 2*Vin)/(2*Vopk) ; 
dn = Dn;

Vapk = Vin/(1-Dn) ;
Vs  = (2*Dn-1)/(Dn*(1-Dn));

Io = Po/Vopk ;
Ila = Io/(1-Dn) ;

Vrms = Vin/sqrt(2);
Iorms = Io /sqrt(2);
Vorms = Vopk/sqrt(2);

Iop = Iorms * sqrt(2) ;
Vop = Vorms * sqrt(2) ; 

Rload = Vopk/Io;   % load resistance (Ohm)
R = Rload ; Rn = R;

Iop = 2*Po/Vop ; 
Vop = 2*Po/Iop ; 
Ilpk = Iop * (1-dn) ;

delta_Vc = Vapk *  delta_V_C_p ; 
delta_Vo = 5/100 ; % 5% 
delta_V_C2 = Vopk * delta_Vo; 
delta_V_C1 = Vin * delta_Vo; 

% dimensions 
Ila_val = Io/(1-Dn);

C1_val = (Io * Dn) / (delta_Vc * fs); Ca = C1_val;     % Capacitor C1 (F)
C2_val = (Io * Dn) / (delta_V_C2 * fs); Cb = C2_val;    % Capacitor C2 (F)

C = Ca - Cb ; 

Iopk = C * Dn / (fs* delta_Vc)  ; 

Iin = Po/Vin ;
IL1_val = Ilpk; I_l1 = IL1_val; I_l2 = I_l1;
Vrms = sqrt(Vin);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;

L1_val = (Vin*(dn))/(delta_I_l1*fs) ; La = L1_val;
L2_val = (Vin*Dn)/(delta_I_l2*fs) ; Lb = L2_val;


% X == [ILaavg; ILbavg; VCaavg; VCbavg]


