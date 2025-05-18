clc
clear all
close all 

% CONVERSOR CUK 
% CHAVE fechada 
% vL1 =Vg
% iC1 =I2
% vL2 =−VC1−Vo 
% iC2 =I2−V2/R
% 
% chave aberta
% vL1 =Vg−VC1
% vL2 =−v2
% iC1 =iL1 
% iC2 =i2−Vo/R



%% REAL VALUES 
Vi_val = 40  ; % Input voltage [V]
Vi=Vi_val;
Vo_val = 80 ; % Output voltage [Vo]
Vo=Vo_val;



D_val = Vo_val / (Vi_val + Vo_val); % SEPIC duty cycle formula


D=D_val;d=D;
fs = 20e3; % [kHz]  1/Ts 
Ts = (1/fs) ; % secs
Vs  = Vi_val*D_val/(1-D_val);
Po = 100 ; % W

%delta_Q = C * delta_Vo
Vrms = sqrt(Vi);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;
Iin = Po/Vi_val ;
IL1_val = Iin; I_l1 = IL1_val;
Io = Po/Vo ;
IL2_val = Io;I_l2= IL2_val;
Rload = Vo/Io;       % load resistance (Ohm)
Ro = Rload ;
delta_I_L_p = 0.2 ; % 20 % percentage
delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;
L2 = (Vi*D)/(delta_I_l2*fs) ;
L1 = (Vi*D)/(delta_I_l1*fs) ;
delta_V_C1_p = 20/100 ; % 20 percentage

delta_V_C1 = Vi * delta_V_C1_p; 
delta_Vo = 5/100 ; % 5% 

delta_V_C2 = Vo * delta_Vo; 
C1_val = (Io * D) / (delta_V_C1 * fs);       % Capacitor C1 (F)

C2_val = (Io * D) / (delta_V_C2 * fs);       % Capacitor C1 (F)
C1=C1_val;C2=C2_val;

delta_Q = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;
%%
num = [
    (C1*L1*Ro*Vi*d - C1*L1*Ro*Vi),
    (L1*Vi*d^2),
    (2*Ro*Vi*d^4 - 4*Ro*Vi*d^3 + Ro*Vi*d^2 + 2*Ro*Vi*d - Ro*Vi)
];

den = [
    (2*C1*C2*L1*L2*Ro*d^3 - 5*C1*C2*L1*L2*Ro*d^2 + 4*C1*C2*L1*L2*Ro*d - C1*C2*L1*L2*Ro),
    (2*C1*L1*L2*d^3 - 5*C1*L1*L2*d^2 + 4*C1*L1*L2*d - C1*L1*L2),
    (C1*L1*Ro - C2*L2*Ro + 13*C1*L1*Ro*d^2 - 12*C1*L1*Ro*d^3 + C2*L1*Ro*d^2 + 4*C1*L1*Ro*d^4 - 4*C2*L1*Ro*d^3 - 14*C2*L2*Ro*d^2 + 5*C2*L1*Ro*d^4 + 16*C2*L2*Ro*d^3 - 2*C2*L1*Ro*d^5 - 9*C2*L2*Ro*d^4 + 2*C2*L2*Ro*d^5 - 6*C1*L1*Ro*d + 6*C2*L2*Ro*d),
    (6*L2*d - L2 + L1*d^2 - 4*L1*d^3 - 14*L2*d^2 + 5*L1*d^4 + 16*L2*d^3 - 2*L1*d^5 - 9*L2*d^4 + 2*L2*d^5),
    (4*Ro*d^6 - 20*Ro*d^5 + 41*Ro*d^4 - 44*Ro*d^3 + 26*Ro*d^2 - 8*Ro*d + Ro)
];

syms s

sys = tf(num, den);
bode(sys)
%%

syms s
Gvo = [((L1*L2*Vi*C1^2*d^2 + 2*L1*L2*Vi*C1*d^2 - L1*L2*Vi*C1*d)*s^3 + (C1^2*L1*Ro*Vi*d - C1*L2*Ro*Vi - 5*C1*L1*Ro*Vi*d^2 - C1*L1*Ro*Vi + 2*C1*L1*Ro*Vi*d^3 - 5*C1*L2*Ro*Vi*d^2 + C1^2*L2*Ro*Vi*d + 2*C1*L2*Ro*Vi*d^3 - 2*C1^2*L1*Ro*Vi*d^2 + C1^2*L1*Ro*Vi*d^3 - 2*C1^2*L2*Ro*Vi*d^2 + C1^2*L2*Ro*Vi*d^3 + 4*C1*L1*Ro*Vi*d + 4*C1*L2*Ro*Vi*d)*s^2 + (2*L2*Vi*d - L1*Vi*d^2 - 6*L2*Vi*d^2 + 2*L1*Vi*d^4 + 6*L2*Vi*d^3 - 2*L2*Vi*d^4 - 2*C1*L2*Vi*d^2 + C1*L1*Vi*d^4 + 5*C1*L2*Vi*d^3 - 3*C1*L2*Vi*d^4 + C1^2*L2*Vi*d^3 - C1^2*L2*Vi*d^4)*s + Ro*Vi*C1^2*d^4 - 2*Ro*Vi*C1^2*d^3 + Ro*Vi*C1^2*d^2 + 4*Ro*Vi*C1*d^4 - 10*Ro*Vi*C1*d^3 + 8*Ro*Vi*C1*d^2 - 2*Ro*Vi*C1*d + 4*Ro*Vi*d^4 - 12*Ro*Vi*d^3 + 13*Ro*Vi*d^2 - 6*Ro*Vi*d + Ro*Vi)/((C2*L1*L2*Ro*C1^2*d^3 - 2*C2*L1*L2*Ro*C1^2*d^2 + C2*L1*L2*Ro*C1^2*d + 2*C2*L1*L2*Ro*C1*d^3 - 5*C2*L1*L2*Ro*C1*d^2 + 4*C2*L1*L2*Ro*C1*d - C2*L1*L2*Ro*C1)*s^4 + (L1*L2*C1^2*d^3 - 2*L1*L2*C1^2*d^2 + L1*L2*C1^2*d + 2*L1*L2*C1*d^3 - 5*L1*L2*C1*d^2 + 4*L1*L2*C1*d - L1*L2*C1)*s^3 + (C1*L1*Ro + C2*L2*Ro + 9*C1*L1*Ro*d^2 - C1^2*L1*Ro*d - 7*C1*L1*Ro*d^3 - C2*L1*Ro*d^2 + 2*C1*L1*Ro*d^4 + 4*C2*L1*Ro*d^3 + 14*C2*L2*Ro*d^2 - 5*C2*L1*Ro*d^4 - 16*C2*L2*Ro*d^3 + 2*C2*L1*Ro*d^5 + 9*C2*L2*Ro*d^4 - 2*C2*L2*Ro*d^5 + 3*C1^2*L1*Ro*d^2 - 3*C1^2*L1*Ro*d^3 + C1^2*L1*Ro*d^4 - 5*C1*L1*Ro*d - 6*C2*L2*Ro*d + C1*C2*L1*Ro*d^3 + 9*C1*C2*L2*Ro*d^2 - 2*C1*C2*L1*Ro*d^4 - 15*C1*C2*L2*Ro*d^3 + C1*C2*L1*Ro*d^5 + 11*C1*C2*L2*Ro*d^4 - 3*C1*C2*L2*Ro*d^5 + C1^2*C2*L2*Ro*d^2 - 3*C1^2*C2*L2*Ro*d^3 + 3*C1^2*C2*L2*Ro*d^4 - C1^2*C2*L2*Ro*d^5 - 2*C1*C2*L2*Ro*d)*s^2 + (L2 - 6*L2*d - L1*d^2 + 4*L1*d^3 + 14*L2*d^2 - 5*L1*d^4 - 16*L2*d^3 + 2*L1*d^5 + 9*L2*d^4 - 2*L2*d^5 + C1^2*L2*d^2 - 3*C1^2*L2*d^3 + 3*C1^2*L2*d^4 - C1^2*L2*d^5 - 2*C1*L2*d + C1*L1*d^3 + 9*C1*L2*d^2 - 2*C1*L1*d^4 - 15*C1*L2*d^3 + C1*L1*d^5 + 11*C1*L2*d^4 - 3*C1*L2*d^5)*s - Ro*C1^2*d^6 + 4*Ro*C1^2*d^5 - 6*Ro*C1^2*d^4 + 4*Ro*C1^2*d^3 - Ro*C1^2*d^2 - 4*Ro*C1*d^6 + 18*Ro*C1*d^5 - 32*Ro*C1*d^4 + 28*Ro*C1*d^3 - 12*Ro*C1*d^2 + 2*Ro*C1*d - 4*Ro*d^6 + 20*Ro*d^5 - 41*Ro*d^4 + 44*Ro*d^3 - 26*Ro*d^2 + 8*Ro*d - Ro), 0, 0, 0]


% Calcular os coeficientes do numerador e denominador
num = [(L1*L2*Vi*C1^2*d^2 + 2*L1*L2*Vi*C1*d^2 - L1*L2*Vi*C1*d), ...
       (C1^2*L1*Ro*Vi*d - C1*L2*Ro*Vi - 5*C1*L1*Ro*Vi*d^2 - C1*L1*Ro*Vi + 2*C1*L1*Ro*Vi*d^3 - 5*C1*L2*Ro*Vi*d^2 + C1^2*L2*Ro*Vi*d + 2*C1*L2*Ro*Vi*d^3 - 2*C1^2*L1*Ro*Vi*d^2 + C1^2*L1*Ro*Vi*d^3 - 2*C1^2*L2*Ro*Vi*d^2 + C1^2*L2*Ro*Vi*d^3 + 4*C1*L1*Ro*Vi*d + 4*C1*L2*Ro*Vi*d), ...
       (2*L2*Vi*d - L1*Vi*d^2 - 6*L2*Vi*d^2 + 2*L1*Vi*d^4 + 6*L2*Vi*d^3 - 2*L2*Vi*d^4 - 2*C1*L2*Vi*d^2 + C1*L1*Vi*d^4 + 5*C1*L2*Vi*d^3 - 3*C1*L2*Vi*d^4 + C1^2*L2*Vi*d^3 - C1^2*L2*Vi*d^4), ...
       (Ro*Vi*C1^2*d^4 - 2*Ro*Vi*C1^2*d^3 + Ro*Vi*C1^2*d^2 + 4*Ro*Vi*C1*d^4 - 10*Ro*Vi*C1*d^3 + 8*Ro*Vi*C1*d^2 - 2*Ro*Vi*C1*d + 4*Ro*Vi*d^4 - 12*Ro*Vi*d^3 + 13*Ro*Vi*d^2 - 6*Ro*Vi*d + Ro*Vi)];

den = [(C2*L1*L2*Ro*C1^2*d^3 - 2*C2*L1*L2*Ro*C1^2*d^2 + C2*L1*L2*Ro*C1^2*d + 2*C2*L1*L2*Ro*C1*d^3 - 5*C2*L1*L2*Ro*C1*d^2 + 4*C2*L1*L2*Ro*C1*d - C2*L1*L2*Ro*C1), ...
       (L1*L2*C1^2*d^3 - 2*L1*L2*C1^2*d^2 + L1*L2*C1^2*d + 2*L1*L2*C1*d^3 - 5*L1*L2*C1*d^2 + 4*L1*L2*C1*d - L1*L2*C1), ...
       (C1*L1*Ro + C2*L2*Ro + 9*C1*L1*Ro*d^2 - C1^2*L1*Ro*d - 7*C1*L1*Ro*d^3 - C2*L1*Ro*d^2 + 2*C1*L1*Ro*d^4 + 4*C2*L1*Ro*d^3 + 14*C2*L2*Ro*d^2 - 5*C2*L1*Ro*d^4 - 16*C2*L2*Ro*d^3 + 2*C2*L1*Ro*d^5 + 9*C2*L2*Ro*d^4 - 2*C2*L2*Ro*d^5 + 3*C1^2*L1*Ro*d^2 - 3*C1^2*L1*Ro*d^3 + C1^2*L1*Ro*d^4 - 5*C1*L1*Ro*d - 6*C2*L2*Ro*d + C1*C2*L1*Ro*d^3 + 9*C1*C2*L2*Ro*d^2 - 2*C1*C2*L1*Ro*d^4 - 15*C1*C2*L2*Ro*d^3 + C1*C2*L1*Ro*d^5 + 11*C1*C2*L2*Ro*d^4 - 3*C1*C2*L2*Ro*d^5 + C1^2*C2*L2*Ro*d^2 - 3*C1^2*C2*L2*Ro*d^3 + 3*C1^2*C2*L2*Ro*d^4 - C1^2*C2*L2*Ro*d^5 - 2*C1*C2*L2*Ro*d), ...
       (L2 - 6*L2*d - L1*d^2 + 4*L1*d^3 + 14*L2*d^2 - 5*L1*d^4 - 16*L2*d^3 + 2*L1*d^5 + 9*L2*d^4 - 2*L2*d^5 + C1^2*L2*d^2 - 3*C1^2*L2*d^3 + 3*C1^2*L2*d^4 - C1^2*L2*d^5 - 2*C1*L2*d + C1*L1*d^3 + 9*C1*L2*d^2 - 2*C1*L1*d^4 - 15*C1*L2*d^3 + C1*L1*d^5 + 11*C1*L2*d^4 - 3*C1*L2*d^5), ...
       (-Ro*C1^2*d^6 + 4*Ro*C1^2*d^5 - 6*Ro*C1^2*d^4 + 4*Ro*C1^2*d^3 - Ro*C1^2*d^2 - 4*Ro*C1*d^6 + 18*Ro*C1*d^5 - 32*Ro*C1*d^4 + 28*Ro*C1*d^3 - 12*Ro*C1*d^2 + 2*Ro*C1*d - 4*Ro*d^6 + 20*Ro*d^5 - 41*Ro*d^4 + 44*Ro*d^3 - 26*Ro*d^2 + 8*Ro*d - Ro)];

% Criar a função de transferência
sys = tf(num, den);
bode(sys)