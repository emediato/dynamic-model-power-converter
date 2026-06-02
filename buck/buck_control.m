% ParÃ¢metros do conversor Buck
Vin = 36;  
Vout = 28;  

C_load_PSU_buck = 10e-6;
L = 470e-6;
R = 100 ;

C = C_load_PSU_buck;
k_v = 1;
V_trip = 1;
k_pwm = 1;

D = Vout / Vin; 

% FunÃ§Ã£o de transferencia de tensÃ£o
s = tf('s');
Gvd = Vin / (1 + (s * L / R) + ((s^2) * L * C));

% FunÃ§Ã£o de transferencia de corrente 
Gid = ((Vin * C * s) + (Vin / R)) / (1 + (s * L / R) + ((s^2) * L * C));
%% Controle analÃ³gico
% Diagrama de Bode da planta de tensÃ£o
figure;
bode(k_pwm*k_v*Gvd);
title('Diagrama de Bode - FunÃ§Ã£o de Transferencia de Tensao');
grid on;

% Diagrama de Bode da planta de corrente
figure;
bode(k_pwm*k_v*Gid);
title('Diagrama de Bode - FunÃ§Ã£o de Transferencia de Corrente');
grid on;

% EspecificaÃ§Ãµes de projeto do controlador PI
Margem_Fase = 50; 
Frequencia_Cruzamento = 100; %%Hz
Frequencia_Cruzamento_rad = 2*pi*Frequencia_Cruzamento;
% CÃ¡lculo do controlador PI 
[mag, phase] = bode(k_pwm*k_v*Gid, Frequencia_Cruzamento_rad);

phase = phase(1) 
mag = mag(1); % Extrai o valor da magnitude

Angulo_PI = Margem_Fase - 180 - phase; 
w_z = Frequencia_Cruzamento_rad / tan(deg2rad(Margem_Fase - 90 - phase));
k_c = Frequencia_Cruzamento_rad / (mag * sqrt(w_z^2 + Frequencia_Cruzamento_rad^2));

PI_Controller = k_c * (1 + s / w_z);

display(PI_Controller);

Constante_de_tempo = 1/(w_z);


%% Controle digital

% ParÃªmetros de entrada do conversor Buck
% clear all;
% close all;
D = 0.5;
Dn = 1-D;
Ls = 240E-6;
Cs = 33E-6;
Ro = 1.92;
Vg = 48;

kv=1;
kpwm=1;
fs=20000;
Ts=1/20e3;
wfiltro=2*pi*2000
MFdeg = 60
Fc = 500 
%Modelo por espaÃ§o de estados
K = [Ls 0;0 Cs];
U = Vg;
A1 = [0 -1; 1 -1/Ro];
A2 = [0 -1; 1 -1/Ro];
B1 = [1; 0];
B2 = [0; 0];
C1 = [1 0; 0 1]; 
C2 = [1 0; 0 1];
E1 = [0; 0];
E2 = [0; 0];

A = D*A1 + Dn*A2;
B = D*B1 + Dn*B2;
C = D*C1 + Dn*C2;
E = D*E1 + Dn*E2;

X = -(inv(A))*B*U;
F = (A1-A2)*X+(B1-B2)*U;
Ap = inv(K)*A;
Bp = inv(K)*B;
Fp = inv(K)*F;

G = (C1-C2)*X+(E1-E2)*Vg;
s = tf([1 0],[1]);

Gy_u = (C*(inv((s*eye(2)-Ap))))*Bp+E;
Gy_d = (C*(inv((s*eye(2)-Ap))))*Fp+G;

Gy1_u = Gy_u(1,:); %Retira a primeira linha 
Gy2_u = Gy_u(2,:);

Gy1_d = Gy_d(1,:)
Gy2_d = Gy_d(2,:)

% FunÃ§Ã£o de transfÃªncia de laÃ§o aberto nÃ£o compensado
Hs=wfiltro/(s+wfiltro)
%Hs=1
FTLAnc=kv*kpwm*Gy2_d*Hs;

%Plot do diagrama de Bode do sistema nÃ£o compensado
%w=10:1:1000000;

% EspecificaÃ§Ãµes para projeto do controlador
wc=2*pi*Fc; %Banda passante do sistema
Mf=MFdeg*pi/180; % Margem de fase projetada

%CÃ¡lculo da mÃ³dulo e da fase na frequÃªncia de interesse
[mod,faseDeg] = bode(FTLAnc,wc);
fase=faseDeg*pi/180;

% CÃ¡lculo dos parÃ¢metros do controlador PI
wz=wc/(tan(Mf-(pi/2)-fase));
kc=wc/(sqrt(wz^2+wc^2)*mod);

% Estrutura do Compensador do tipo PI
Cv=kc*tf([1 wz],[1 0]);

%FunÃ§Ã£o de transferÃªncia de laÃ§o aberto compensado
FTLAc=FTLAnc*Cv;

%Resposta ao degrau do sistema compensado
step(FTLAc/(FTLAc*Hs+1));

%DiscretizaÃ§Ã£o por ZOH da planta
Gz2=c2d(FTLAnc*Hs,Ts);
z=tf('z',Ts);
Gz=Gz2/z;
%DiscretizaÃ§Ã£o por ZOH do compensador
%Somente para fins de visualizaÃ§Ã£o
Cz=c2d(Cv,Ts);
hold on
% Resposta ao degrau do sistema discretizado por ZOH
% Importante notar que as respostas divergem
step(Cz*Gz/(Cz*Gz+1));

% ConversÃ£o do plano Z para o plano W
Gw = d2c(Gz,'tustin');

%Pre-warping
wc2 = (2/Ts)*tan(wc*Ts/2);

%MÃ³dulo e fase do bode discreto
[mag2,phase2] = bode(Gw,wc2);
%Se phase 2 for muito do que zero deve-se diminuir de 2*pi
phaseRad2 = phase2*pi/180 - 2*pi;

%Calculo dos parÃ¢metros do controlador PI
wz2=wc2/(tan(Mf-(pi/2)-phaseRad2));
kc2=wc2/(sqrt(wz2^2+wc2^2)*mag2);

% Estrutura de compensador em W
Cv2=kc2*tf([1 wz2],[1 0]);

%DiscretizaÃ§Ã£o do compensador para Z
Cz2 = c2d(Cv2,Ts,'tustin');

% Resposta do compensador discretizado
% Essa resposta deve ser muito prÃ³xima a resposta analÃ³gica
step(Cz2*Gz/(Cz2*Gz+1));
