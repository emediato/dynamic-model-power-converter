%% ==============================
% Projeto de Controlador de Corrente - Conversor Buck
% Comparação entre planta simplificada e planta realista
% ==============================

clear; close all; clc;

%% 1. Parâmetros do conversor
Vin = 100;      % V
Ro = 6;         % Ohm
Lo = 2e-3;      % H
Co = 10e-6;     % F
Io = Vin/Ro;    % A (valor de referência)

%% 2. Ganhos do sistema
ki = 0.33;      % V/A (sensor de corrente)
kpwm = 1;       % modulador PWM (ganho unitário)

%% 3. Funções de transferência da planta (corrente)
s = tf('s');

% Planta A: integrador puro
GA = Vin / (s * Lo);

% Planta B: modelo completo com carga e capacitor
% Gid(s) = Vo(s)/d(s) / Ro? Na verdade, modelo corrente:
% Para buck: Gid(s) = Vin/(s*Lo) * (1/(1 + s*Ro*Co))?
% Vamos usar modelo mais realista:
GB = Vin / (s*Lo + Ro/(1 + s*Ro*Co));

%% 4. Função de transferência de laço aberto não compensada (FTLAnc)
FTLAnc_A = GA * ki * kpwm;
FTLAnc_B = GB * ki * kpwm;

%% 5. Especificações de projeto
fci = 2000;             % Hz
wci = 2*pi*fci;         % rad/s
Mphi = 60 * pi/180;     % 60° margem de fase

%% 6. Projeto do compensador (avanço de fase) - para planta A
[magA, phaseA] = bode(FTLAnc_A, wci);
phaseA_deg = phaseA;

phi_lead = rad2deg(Mphi) - 90 - phaseA_deg;
phi_lead_rad = deg2rad(phi_lead);

wzA = wci / tan(phi_lead_rad + pi/2);
kcA = wci / (magA * sqrt(wci^2 + wzA^2));

CiA = kcA * (s + wzA) / s;

%% 7. Projeto do compensador para planta B
[magB, phaseB] = bode(FTLAnc_B, wci);
phaseB_deg = phaseB;

phi_leadB = rad2deg(Mphi) - 90 - phaseB_deg;
phi_lead_radB = deg2rad(phi_leadB);

wzB = wci / tan(phi_lead_radB + pi/2);
kcB = wci / (magB * sqrt(wci^2 + wzB^2));

CiB = kcB * (s + wzB) / s;

%% 8. FTLA compensada
FTLAc_A = FTLAnc_A * CiA;
FTLAc_B = FTLAnc_B * CiB;

%% 9. Margens
[Gm_A, Pm_A, ~, ~] = margin(FTLAc_A);
[Gm_B, Pm_B, ~, ~] = margin(FTLAc_B);

%% 10. Malha fechada
MF_A = feedback(FTLAc_A, 1);
MF_B = feedback(FTLAc_B, 1);

%% 11. Resposta ao degrau
t = 0:1e-5:2e-3;  % 2ms
[yA, tA] = step(MF_A, t);
[yB, tB] = step(MF_B, t);

%% 12. Resultados
fprintf('========== PLANTA A (INTEGRADOR PURO) ==========\n');
fprintf('Ganho kc: %.3f\n', kcA);
fprintf('Zero wz: %.2f rad/s\n', wzA);
fprintf('Margem de fase obtida: %.2f graus\n', Pm_A);
fprintf('Margem de ganho obtida: %.2f dB\n', 20*log10(Gm_A));
fprintf('================================================\n\n');

fprintf('========== PLANTA B (MODELO COMPLETO) ==========\n');
fprintf('Ganho kc: %.3f\n', kcB);
fprintf('Zero wz: %.2f rad/s\n', wzB);
fprintf('Margem de fase obtida: %.2f graus\n', Pm_B);
fprintf('Margem de ganho obtida: %.2f dB\n', 20*log10(Gm_B));
fprintf('================================================\n\n');

%% 13. Gráficos comparativos
figure;
subplot(2,1,1);
bode(FTLAnc_A, FTLAc_A, {2*pi*100, 2*pi*10000});
legend('FTLA não compensada (A)', 'FTLA compensada (A)', 'Location', 'southwest');
title('Bode - Planta A (integradora)');
grid on;

subplot(2,1,2);
bode(FTLAnc_B, FTLAc_B, {2*pi*100, 2*pi*10000});
legend('FTLA não compensada (B)', 'FTLA compensada (B)', 'Location', 'southwest');
title('Bode - Planta B (completa)');
grid on;

figure;
subplot(2,1,1);
step(MF_A, t);
title('Resposta ao degrau - Planta A');
xlabel('Tempo (s)'); ylabel('Corrente (A)');
grid on;

subplot(2,1,2);
step(MF_B, t);
title('Resposta ao degrau - Planta B');
xlabel('Tempo (s)'); ylabel('Corrente (A)');
grid on;

%% 14. Comparação direta das FTLA compensadas
figure;
bode(FTLAc_A, FTLAc_B, {2*pi*100, 2*pi*10000});
legend('FTLA compensada - Planta A', 'FTLA compensada - Planta B');
title('Comparação das FTLA compensadas');
grid on;
