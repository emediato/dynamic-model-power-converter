
clc
clear
close all


%% Variáveis de estado
syms d Vin C1 C2 L1 L2 Ro s

U = [Vin; 0; 0; 0];

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

% CC MODEL - operation point behaviour
X = - simplify(inv(AA)*BB)*U ;

% Transfer function
E1 = ( AA1 - AA2)*X;
E2 = simplify((BB1-BB2)*U);
E = simplify(E1+E2);
CC1 = [0 0 0 1];

%colect isola variavel, ex em funcao de s
Gvo = collect(simplify(CC1*inv(s*eye(4)-AA)*E), s)

% Ca = [0 0 0 1];
% Gvc1 = collect(simplify(Ca*inv(s*eye(4)-AA)*BB,s))
CC1 = [0 0 0 1]; % Vo/d
CC2 = [1 0 0 0]; % i_L/d


% y / d
G_vo_d = CC1 / (s * eye(4) - AA) * [E];
G_vo_d = collect(simplify(G_vo_d), s)

% IL / d
G_iL_d = CC2 / (s * eye(4) - AA) * [E];
G_iL_d = collect(simplify(G_iL_d), s)

%% REAL VALUES  Conversor
Vi = 48;           % Tensão de entrada (V)
Vo = 72;            % Tensão de saída (V)
delta_iL = 0.2; delta_I_L_p = delta_iL;    % Variação da corrente (relativa)
delta_vC1 = 0.2; delta_V_C1_p = delta_vC1;   % Variação da tensão C1 (relativa)
delta_Vo = 0.05;    % Variação da tensão de saída (relativa)
fs = 50e3;          % Frequência de comutação (Hz)
Po = 300;           % Potência de saída (W)

D_val = Vo / (Vi + Vo); % Cuk/sepic duty cycle formula
D=D_val;

Ts = (1/fs) ; % secs
Vs  = Vi*D_val/(1-D_val);

Io = Po/Vo ;
Iin = Po/Vi ;
IL1_val = Iin; I_l1 = IL1_val;

IL2_val = Io;I_l2= IL2_val;
Rload = Vo/Io;       % load resistance (Ohm)
Rload2 = (Vo^2) / Po;

%delta_Q = C * delta_Vo
Vrms = sqrt(Vi);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;

L1_val = (Vi*(D))/(delta_I_l1*fs) ;
L2_val = (Vo*(1-D))/(delta_I_l2*fs) ;

delta_V_C1 = Vi * delta_V_C1_p; 
delta_V_C2 = Vo * delta_Vo; 

C1_val = (Io * D) / (delta_V_C1 * fs);       % Capacitor C1 (F)
C2_val = (Vo * (1-D)) / (8 * delta_V_C2 * fs^2 * L2_val)       % Capacitor C2 (F)


delta_Q = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;

%%
% Definindo os valores reais
valores = [C1, C2, L1, L2, Ro, Vin, d]; %simbolic
numeros = [C1_val, C2_val, L1_val, L2_val,  Rload, Vi,  D]; % Exemplo de valores

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
