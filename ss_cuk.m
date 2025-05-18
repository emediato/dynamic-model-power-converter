
clc
clear
close all

% 1a e 2a etapa
% ic = il - io
%filtrar ondulacao da corrente do indutor

% il - io / delta_il/2 
% area do triangulo : b*2/2 .
% deltaq = c * DELTA_VO
% DELTAQ = (TS/2 * DELTA_IL/2)/2

% isolar C

%% Variáveis de estado
syms d Vi C1 C2 L1 L2 Ro s

U = [Vi 0 0 0];

AA1 = [0 0 0 0 ;
    0 0 1/L2 -1/L2; 
    0 -1/C1 0 0;
    0 1/C2 0 -1/(Ro*C2)];

AA2 = [0 0 -1/L1 0;
        0 0 0 -1/L2;
        1/C1 0 0 0;
        0 1/C2 0 -1/(Ro*C2)];

AA = simplify(AA1*d+AA2*(1-d));
A = AA;


BB1= [1/L1  0     0      0;
    0      0     0      0;
    0      0     0      0;
    0      0     0      0];
BB2=BB1;

BB = simplify((BB1*d)+(BB2*(1-d)));
B=BB;

%% CC MODEL - operation point behaviour
X = - simplify(inv(AA)*BB)*U ;

% Transfer function
E1 = ( AA1 - AA2)*X;

E2 = simplify((BB1-BB2)*U);

E = simplify(E1+E2);

C1 = [0 0 0 1];

%colect isola variavel, ex em funcao de s
Gvo = collect(simplify(C1*inv(s*eye(4)-AA)*E), s)

% GVo = [((L1*L2*Vi*C1^2*d^2 + 2*L1*L2*Vi*C1*d^2 - L1*L2*Vi*C1*d)*s^3 + (C1^2*L1*Ro*Vi*d - C1*L2*Ro*Vi - 5*C1*L1*Ro*Vi*d^2 - C1*L1*Ro*Vi + 2*C1*L1*Ro*Vi*d^3 - 5*C1*L2*Ro*Vi*d^2 + C1^2*L2*Ro*Vi*d + 2*C1*L2*Ro*Vi*d^3 - 2*C1^2*L1*Ro*Vi*d^2 + C1^2*L1*Ro*Vi*d^3 - 2*C1^2*L2*Ro*Vi*d^2 + C1^2*L2*Ro*Vi*d^3 + 4*C1*L1*Ro*Vi*d + 4*C1*L2*Ro*Vi*d)*s^2 + (2*L2*Vi*d - L1*Vi*d^2 - 6*L2*Vi*d^2 + 2*L1*Vi*d^4 + 6*L2*Vi*d^3 - 2*L2*Vi*d^4 - 2*C1*L2*Vi*d^2 + C1*L1*Vi*d^4 + 5*C1*L2*Vi*d^3 - 3*C1*L2*Vi*d^4 + C1^2*L2*Vi*d^3 - C1^2*L2*Vi*d^4)*s + Ro*Vi*C1^2*d^4 - 2*Ro*Vi*C1^2*d^3 + Ro*Vi*C1^2*d^2 + 4*Ro*Vi*C1*d^4 - 10*Ro*Vi*C1*d^3 + 8*Ro*Vi*C1*d^2 - 2*Ro*Vi*C1*d + 4*Ro*Vi*d^4 - 12*Ro*Vi*d^3 + 13*Ro*Vi*d^2 - 6*Ro*Vi*d + Ro*Vi)/((C2*L1*L2*Ro*C1^2*d^3 - 2*C2*L1*L2*Ro*C1^2*d^2 + C2*L1*L2*Ro*C1^2*d + 2*C2*L1*L2*Ro*C1*d^3 - 5*C2*L1*L2*Ro*C1*d^2 + 4*C2*L1*L2*Ro*C1*d - C2*L1*L2*Ro*C1)*s^4 + (L1*L2*C1^2*d^3 - 2*L1*L2*C1^2*d^2 + L1*L2*C1^2*d + 2*L1*L2*C1*d^3 - 5*L1*L2*C1*d^2 + 4*L1*L2*C1*d - L1*L2*C1)*s^3 + (C1*L1*Ro + C2*L2*Ro + 9*C1*L1*Ro*d^2 - C1^2*L1*Ro*d - 7*C1*L1*Ro*d^3 - C2*L1*Ro*d^2 + 2*C1*L1*Ro*d^4 + 4*C2*L1*Ro*d^3 + 14*C2*L2*Ro*d^2 - 5*C2*L1*Ro*d^4 - 16*C2*L2*Ro*d^3 + 2*C2*L1*Ro*d^5 + 9*C2*L2*Ro*d^4 - 2*C2*L2*Ro*d^5 + 3*C1^2*L1*Ro*d^2 - 3*C1^2*L1*Ro*d^3 + C1^2*L1*Ro*d^4 - 5*C1*L1*Ro*d - 6*C2*L2*Ro*d + C1*C2*L1*Ro*d^3 + 9*C1*C2*L2*Ro*d^2 - 2*C1*C2*L1*Ro*d^4 - 15*C1*C2*L2*Ro*d^3 + C1*C2*L1*Ro*d^5 + 11*C1*C2*L2*Ro*d^4 - 3*C1*C2*L2*Ro*d^5 + C1^2*C2*L2*Ro*d^2 - 3*C1^2*C2*L2*Ro*d^3 + 3*C1^2*C2*L2*Ro*d^4 - C1^2*C2*L2*Ro*d^5 - 2*C1*C2*L2*Ro*d)*s^2 + (L2 - 6*L2*d - L1*d^2 + 4*L1*d^3 + 14*L2*d^2 - 5*L1*d^4 - 16*L2*d^3 + 2*L1*d^5 + 9*L2*d^4 - 2*L2*d^5 + C1^2*L2*d^2 - 3*C1^2*L2*d^3 + 3*C1^2*L2*d^4 - C1^2*L2*d^5 - 2*C1*L2*d + C1*L1*d^3 + 9*C1*L2*d^2 - 2*C1*L1*d^4 - 15*C1*L2*d^3 + C1*L1*d^5 + 11*C1*L2*d^4 - 3*C1*L2*d^5)*s - Ro*C1^2*d^6 + 4*Ro*C1^2*d^5 - 6*Ro*C1^2*d^4 + 4*Ro*C1^2*d^3 - Ro*C1^2*d^2 - 4*Ro*C1*d^6 + 18*Ro*C1*d^5 - 32*Ro*C1*d^4 + 28*Ro*C1*d^3 - 12*Ro*C1*d^2 + 2*Ro*C1*d - 4*Ro*d^6 + 20*Ro*d^5 - 41*Ro*d^4 + 44*Ro*d^3 - 26*Ro*d^2 + 8*Ro*d - Ro), 0, 0, 0] ;
%%[((C1*L1*Ro*Vi - C1*L1*Ro*Vi*d)*s^2 + (-L1*Vi*d^2)*s + Ro*Vi*d^2 - 2*Ro*Vi*d + Ro*Vi)/((C1*C2*L1*L2*Ro*d^2 - 2*C1*C2*L1*L2*Ro*d + C1*C2*L1*L2*Ro)*s^4 + (C1*L1*L2*d^2 - 2*C1*L1*L2*d + C1*L1*L2)*s^3 + (C1*L1*Ro + C2*L2*Ro + C1*L1*Ro*d^2 + C2*L1*Ro*d^2 - 2*C2*L1*Ro*d^3 + 6*C2*L2*Ro*d^2 + C2*L1*Ro*d^4 - 4*C2*L2*Ro*d^3 + C2*L2*Ro*d^4 - 2*C1*L1*Ro*d - 4*C2*L2*Ro*d)*s^2 + (L2 - 4*L2*d + L1*d^2 - 2*L1*d^3 + 6*L2*d^2 + L1*d^4 - 4*L2*d^3 + L2*d^4)*s + Ro*d^4 - 4*Ro*d^3 + 6*Ro*d^2 - 4*Ro*d + Ro), 0, 0, 0]

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

Io = Po/Vo ;
Iin = Po/Vi_val ;
IL1_val = Iin; I_l1 = IL1_val;

IL2_val = Io;I_l2= IL2_val;
Rload = Vo/Io;       % load resistance (Ohm)
Ro = Rload ;

%delta_Q = C * delta_Vo
Vrms = sqrt(Vi);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

delta_I_L_p = 0.2 ; % 20 % percentage

delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;

L1 = (Vi*(D))/(delta_I_l1*fs) ;
L2 = (Vo*(1-D))/(delta_I_l2*fs) ;

delta_V_C1_p = 20/100 ; % 20 percentage
delta_V_C1 = Vi * delta_V_C1_p; 
delta_Vo = 5/100 ; % 5% 
delta_V_C2 = Vo * delta_Vo; 

C1_val = (Io * D) / (delta_V_C1 * fs);       % Capacitor C1 (F)
C2_val = (Vo * (1-D)) / (8 * delta_V_C2 * fs^2 * L2)       % Capacitor C2 (F)
C1=C1_val;C2=C2_val;

delta_Q = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;
%%
% Calcular os coeficientes do numerador e denominador% Termo s² (a₂)
%%[((C1*L1*Ro*Vi - C1*L1*Ro*Vi*d)*s^2 + (-L1*Vi*d^2)*s + Ro*Vi*d^2 - 2*Ro*Vi*d + Ro*Vi)/((C1*C2*L1*L2*Ro*d^2 - 2*C1*C2*L1*L2*Ro*d + C1*C2*L1*L2*Ro)*s^4 + (C1*L1*L2*d^2 - 2*C1*L1*L2*d + C1*L1*L2)*s^3 + (C1*L1*Ro + C2*L2*Ro + C1*L1*Ro*d^2 + C2*L1*Ro*d^2 - 2*C2*L1*Ro*d^3 + 6*C2*L2*Ro*d^2 + C2*L1*Ro*d^4 - 4*C2*L2*Ro*d^3 + C2*L2*Ro*d^4 - 2*C1*L1*Ro*d - 4*C2*L2*Ro*d)*s^2 + (L2 - 4*L2*d + L1*d^2 - 2*L1*d^3 + 6*L2*d^2 + L1*d^4 - 4*L2*d^3 + L2*d^4)*s + Ro*d^4 - 4*Ro*d^3 + 6*Ro*d^2 - 4*Ro*d + Ro), 0, 0, 0]

%%[((C1*L1*Ro*Vi - C1*L1*Ro*Vi*d)*s^2 + (-L1*Vi*d^2)*s + Ro*Vi*d^2 - 2*Ro*Vi*d + Ro*Vi)
% --------------------------------
% /((C1*C2*L1*L2*Ro*d^2 - 2*C1*C2*L1*L2*Ro*d + C1*C2*L1*L2*Ro)*s^4 + (C1*L1*L2*d^2 - 2*C1*L1*L2*d + C1*L1*L2)*s^3 + (C1*L1*Ro + C2*L2*Ro + C1*L1*Ro*d^2 + C2*L1*Ro*d^2 - 2*C2*L1*Ro*d^3 + 6*C2*L2*Ro*d^2 + C2*L1*Ro*d^4 - 4*C2*L2*Ro*d^3 + C2*L2*Ro*d^4 - 2*C1*L1*Ro*d - 4*C2*L2*Ro*d)*s^2 + (L2 - 4*L2*d + L1*d^2 - 2*L1*d^3 + 6*L2*d^2 + L1*d^4 - 4*L2*d^3 + L2*d^4)*s + Ro*d^4 - 4*Ro*d^3 + 6*Ro*d^2 - 4*Ro*d + Ro), 0, 0, 0]

%[(s^2/84375 - (32*s)/675 + 2560/9)/((521666913674169*s^4)/2535301200456458802993406410752 + (636800431731163*s^3)/77371252455336267181195264 + (244863447259907*s^2)/2361183241434822606848 + (2*s)/10125 + 64/81), 0, 0, 0]

% Coeficientes do numerador (ordem decrescente: s^2, s^1, s^0)
% (s^2/84375 - (32*s)/675 + 2560/9)

ns2 =  1/84375;
ns1 =  -32/675;
ns0 =  2560/9; 
num = [ns2, ns1, 284.444444];

% Coeficientes do denominador (ordem decrescente: s^4, s^3, s^2, s^1, s^0)
% /((521666913674169*s^4)/2535301200456458802993406410752 + (636800431731163*s^3)/77371252455336267181195264 + (244863447259907*s^2)/2361183241434822606848 + (2*s)/10125 + 64/81), 0, 0, 0]
% 521666913674169/2535301200456458802993406410752 
% + (636800431731163*s^3)/77371252455336267181195264 
% + (244863447259907*s^2)/2361183241434822606848 
% + (2*s)/10125 
% + 64/81), 0, 0, 0]

dens4 = 521666913674169/2535301200456458802993406410752 ;
dens3 = 636800431731163/77371252455336267181195264 ;
dens2 = 244863447259907/2361183241434822606848 ; 
dens1 = 2/10125 ;
dens0= 64/81;
den = [dens4,dens3,dens2,dens1,dens0];

G = tf(num, den);
% Criando a função de transferência
sys = tf(num, den);

% Exibir a função de transferência
disp('Função de transferência:');
disp(sys);
%%
bode(sys)
%%
step(sys)

