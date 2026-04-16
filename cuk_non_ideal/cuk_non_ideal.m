clc
clear
close all
% dec2hex(bin2dec('01010111'))
% dec2bin(hex2dec('0x90')) '10010000'
% Variáveis de estado
syms C1 C2 L1 L2 Ro Vo il1 il2 dil1 dil2 vs Vf Vi rL rc ro rs rd rc1 rc2 rL1 rL2 Vc1 Vc2 Vco dvc1 dvo s d
eqn1 = L1*dil1== Vi -il1*(rL+rs) - rs*il2 ;
eqn2 = L2*dil2 == -Vo +Vc1 -il1*(rs) -il2*(rL+rc+rs) ;
eqn3 = C1*dvc1 == il2 ;
eqn4 = C2*dvo == -Vo/Ro + il2 ;
[s1_dil1, s1_dil2, s1_dvc1, s1_dvo] = solve(eqn1, eqn2, eqn3, eqn4, dil1, dil2, dvc1, dvo);
s1_dil1_isolated = collect(simplify(s1_dil1), [il1 il2 Vc1 Vi] )
s1_dil2_isolated = collect(simplify(s1_dil2), [il1 il2 Vc1 Vf])
s1_dvc1_isolated = collect(simplify(s1_dvc1), Vc1) % // 0
s1_dvo_isolated = collect(simplify(s1_dvo), Vo) % // 0
% 2a etapa
eqn5 = L1*dil1 == Vi - Vf - Vc1 * (-il1*(rL+rc+rd)) - rd*il2 ;
eqn6 = L2*dil2 == - Vo - Vf - il1*(rL+rd+rc) - rd*il2 ;
eqn7 = C1*dvc1 == il1 ;
eqn8 = C2*dvo == il2 - Vo/Ro ;
[s2_dil, s2_dil2, s2_dvc1, s2_dvo] = solve(eqn5, eqn6, eqn7, eqn8, dil1, dil2, dvc1, dvo);
s2_dl1_isolated = collect(simplify(s2_dil), [il1 il2 Vi Vc1])
s2_dl2_isolated = collect(simplify(s2_dil2), [il1 il2 Vo ])
s2_dvc1_isolated = collect(simplify(s2_dvc1), Vc1) % // 0
s2_dvo_isolated = collect(simplify(s2_dvo), [il1 il2 Vo]) % // 0
%%
U = [Vi; 0; Vf; 0] ;
AA1=[-(rs+rL)/L1, -rs/L1, 0, 0;
-(rs)/L2, (-rL+rc+rs)/L2, 1/L2, -1/L2;
0, 1/C1, 0, 0;
0, 1/C2, 0, -1/(C2*Ro)];
AA2=[(rL+rc+rd)/L1, -rd/L1, (rL+rc+rd)/L1, 0;
-rd/L2, (-(Ro + rL + rd))/L2, 1/L2, 0;
1/C1, 0, 0, 0;
0, +1/C2, 0, -1/(C2*Ro)];
AA = simplify(AA1*d+AA2*(1-d));
A = AA;
BB1 = [1/L1 0 0 0 ; 0 0 0 0; 0 0 0 0; 0 0 0 0];
BB2 = [1/L1 -1/L2 0 0 ;0 -1/L2 0 0; 0 0 0 0; 0 0 0 0];
BB = simplify((BB1*d)+(BB2*(1-d)));
BB = simplify((BB1*d)+(BB2*(1-d)));
B=BB;
% CC MODEL - operation point behaviour
X = - simplify(inv(AA)*BB)*U ;
% Transfer function
E1 = ( AA1 - AA2)*X;
E2 = simplify((BB1-BB2)*U);
E = simplify(E1+E2);
CC1 = [0 0 0 1];
%colect isola variavel em funcao de s
Gvo = collect(simplify(CC1*inv(s*eye(4)-AA)*E), s)
%% MODELO CC
% regime permanente
X=-simplify(inv(A)*B)*U ;
% FT
E1 = simplify((AA1-AA2)*X);
E2 = simplify((BB1-BB2)*U);
E = simplify(E1+E2);
%
Ca = [0 1 0 0];
Cb = [0 0 1 0];
G1 = Ca / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);
G2 = Cb / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);
% y / d
G = (simplify(G1+G2));
Gvod = collect(simplify(G), s);
CC1 = [0 0 0 1]; % Vo/d
CC2 = [1 0 0 0]; % i_L/d
% y / d
G_vo_d = CC1 / (s * eye(4) - AA) * [E];
G_vo_d = collect(simplify(G_vo_d), s)
% IL / d
G_iL_d = CC2 / (s * eye(4) - AA) * [E];
G_iL_d = collect(simplify(G_iL_d), s)
%% REAL VALUES Conversor
%% nominal values
Vin = 48; % Tensão de entrada (V)
Vo = 72; % Tensão de saída (V)
delta_iL = 0.2; delta_I_L_p = delta_iL; % Variação da corrente (relativa)
delta_vC1 = 0.2; delta_V_C1_p = delta_vC1; % Variação da tensão C1 (relativa)
delta_Vo = 0.05; % Variação da tensão de saída (relativa)
fs = 50e3; % Frequência de comutação (Hz)
Dn = Vo/(Vin+Vo) ; dn = Dn;
Vs = Vin*Dn/(1-Dn);
Po = 250 ; % W
rLn = 0.1 ; %100 %mΩ
rcn= 0.010 ; % 10 mΩ
rsn = .08 ;% 80 mΩ
rdn = .08 ;% 80 mΩ
VFn = 1 ; % V
D_val = Vo / (Vin + Vo); % Cuk/sepic duty cycle formula
D=D_val;
Ts = (1/fs) ; % secs
Vs = Vin*D_val/(1-D_val);
Io = Po/Vo ;
Iin = Po/Vin ;
IL1_val = Iin; I_l1 = IL1_val;
IL2_val = Io;I_l2= IL2_val;
Rload = Vo/Io; % load resistance (Ohm)
Rload2 = (Vo^2) / Po;
%delta_Q = C * delta_Vo
Vrms = sqrt(Vin);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;
delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;
L1_val = (Vin*(D))/(delta_I_l1*fs) ;
L2_val = (Vo*(1-D))/(delta_I_l2*fs) ;
delta_V_C1 = Vin * delta_V_C1_p;
delta_V_C2 = Vo * delta_Vo;
C1_val = (Io * D) / (delta_V_C1 * fs); % Capacitor C1 (F)
C2_val = (Vo * (1-D)) / (8 * delta_V_C2 * fs^2 * L2_val); % Capacitor C2 (F)
delta_Q = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;
%%
vars = [Ro L1 L2 C1 C2 Vi d rL rc rs rd Vf ];
vals = [Rload L1_val L2_val C1_val C2_val Vin Dn rLn rcn rsn rdn VFn ];
Gn_vod = subs(Gvod, vars, vals);
Gn_il = subs(G_iL_d, vars, vals);
[Gn_num, Gn_den] = numden(Gn_il);
tf_num = double(flip(coeffs(Gn_num)));
tf_den = double(flip(coeffs(Gn_den)));
Gntf_il = tf(tf_num, tf_den);
step(Gntf_il)
%%
%
% Definindo os valores reais
valores = [C1 C2 L1 L2 Ro Vi d rL rc rs rd Vf]; %simbolic
numeros = [C1_val C2_val L1_val L2_val Rload Vin D rLn rcn rsn rdn VFn]; % valores nominais
% Substituição (exceto o 's')
Gvo_num_sym = subs(Gvo, valores, numeros);
G_iL_d_num_sym = subs(G_iL_d, valores, numeros);
% Separa numerador e denominador do resultado simbólico
[Num, Den] = numden(Gvo_num_sym);
[Num2, Den2] = numden(G_iL_d_num_sym);
% Converte os polinômios simbólicos em vetores numéricos de coeficientes
% O MATLAB extrai na ordem: [s^n, s^(n-1), ..., s^0]
num_vetor = double(sym2poly(Num));
den_vetor = double(sym2poly(Den));
num_vetor2 = double(sym2poly(Num2));
den_vetor2 = double(sym2poly(Den2));
% Se o ganho ficar muito alto no numerador, pode ser necessário normalizar
% dividindo ambos pelo primeiro termo do denominador:
num_vetor = num_vetor / den_vetor(1);
den_vetor = den_vetor / den_vetor(1);
num_vetor2 = num_vetor2 / den_vetor2(1);
den_vetor2 = den_vetor2 / den_vetor2(1);
% Cria a Função de Transferência Numérica
sys = tf(num_vetor, den_vetor);
sys_iL = tf(num_vetor2, den_vetor2);
%sys = tf(num, den);
% Exibir a função de transferência
disp('Função de transferência: Gvod');
disp(sys);
% Exibir a função de transferência
disp('Função de transferência: GiLd');
disp(sys_iL);
%%
figure('Name', 'Análise Dinâmica do Conversor Cuk - G_{vd}(s)', 'NumberTitle', 'off');
% Step Response
subplot(2,1,1);
step(sys);
grid on;
title('Resposta ao Degrau Unitário - G_{vd}(s)');
xlabel('Tempo (s)');
% Diagrama de Bode
subplot(2,1,2);
bode(sys);
grid on;
title('Diagrama de Bode - G_{vd}(s)');
% Cálculo de Margens de Estabilidade
[gm, pm, wcg, wcp] = margin(sys);
fprintf('\n--- Análise de Estabilidade ---\n');
fprintf('Margem de Fase: %.2f graus em %.2f rad/s\n', pm, wcp);
fprintf('Margem de Ganho: %.2f dB em %.2f rad/s\n', 20*log10(gm), wcg);
%%
figure('Name', 'Análise Dinâmica do Conversor Cuk G_{iLd}', 'NumberTitle', 'off');
% Step Response
subplot(2,1,1);
step(sys_iL);
grid on;
title('Resposta ao Degrau Unitário - G_{iLd}(s)');
xlabel('Tempo (s)');
% Diagrama de Bode
subplot(2,1,2);
bode(sys_iL);
grid on;
title('Diagrama de Bode - G_{iLd}(s)');
% Cálculo de Margens de Estabilidade
[gm, pm, wcg, wcp] = margin(sys_iL);
fprintf('\n--- Análise de Estabilidade ---\n');
fprintf('Margem de Fase: %.2f graus em %.2f rad/s\n', pm, wcp);
fprintf('Margem de Ganho: %.2f dB em %.2f rad/s\n', 20*log10(gm), wcg);
