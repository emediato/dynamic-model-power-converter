%% ========================================================================
clear; close all; clc;

%% 1. System Parameters
Vin = 100;      % V
Ro = 6;         % Ohm
Lo = 2e-3;      % H
Co = 10e-6;     % F
fs = 40e3;      % Hz

% Ganhos de Sensores e PWM
ki = 0.33;      % V/A (Sensor de Corrente)
kv = 0.055;     % V/V (Sensor de Tensão: 3.3V / 60V)
kpwm = 1;       % Ganho do Modulador

s = tf('s');

%% 2. Inner Loop: Current Control 

Gid = Vin / (s*Lo + Ro/(1 + s*Ro*Co));
% Especificações de Corrente
fci = 2000;             % Crossover 2kHz
wci = 2*pi*fci;
Mphi_i = 60 * pi/180;   % Margem de fase 60 graus

% Sintonia do Compensador Ci(s)
[mag_i, phase_i] = bode(Gid * ki * kpwm, wci);
mag_i = squeeze(mag_i); phase_i = squeeze(phase_i);


% wz_i é o omega do zero do controlador
wzi = wci / tan(Mphi_i - pi/2 - deg2rad(phase_i));
kci = wci / (mag_i * sqrt(wci^2 + wzi^2));
Ci = kci * (s + wzi) / s;

tau_i = 1 / wzi;           % Constante de tempo (s)

% Malha Aberta e Fechada de Corrente
FTLA_i = Gid * ki * kpwm * Ci;

% FTMF_i = feedback(FTLA_i, 1)
FTMF_i = Ci * kpwm * Gid / ( 1 + ki * (Ci * kpwm * Gid)) ;
[mag_ftmf_i, phase_ftmf_i] = bode(FTMF_i, wzi);
mag_ftmf_i = squeeze(mag_ftmf_i); phase_ftmf_i = squeeze(phase_ftmf_i);


%% 3. Outer Loop: Voltage Control (Cascaded)
%%Projeto da Malha Externa (Tensão)
% Planta de Tensão: Gvi(s) = Vo(s) / iL(s)
Gvi = Ro / (s * Ro * Co + 1);

fcv = 200;                      % Frequência de cruzamento desejada (Hz)
omega_cv = 2 * pi * fcv;
wcv = 2*pi*fcv;       % Frequência de cruzamento (rad/s)
Mphi_v = 90 * pi/180;           % Margem de Fase (rad)


FTLA_ncv = Gvi * kv * FTMF_i;   % Malha aberta não compensada
[mag_v, phase_v] = bode(FTLA_ncv, omega_cv);
mag_v = squeeze(mag_v); phase_v = squeeze(phase_v);
 
% Cálculo do zero do controlador (ωz)
stan_arg_v = tan(Mphi_v - pi/2 - phase_v);
wzv = omega_cv / stan_arg_v
kcv = wcv / (mag_v * sqrt(wcv^2 + wzv^2));
tau_v = 1 / wzv


Cv = kcv * (s + wzv) / s;

% Cálculo do Zero e Tau de Tensão
omega_zv = omega_cv / tan(Mphi_v - pi/2 - deg2rad(phase_v));
tau_v = 1 / omega_zv;           % Constante de tempo (s)
kc_v = omega_cv / (mag_v * sqrt(omega_cv^2 + omega_zv^2));

fprintf('--- PARÂMETROS DA MALHA DE CORRENTE ---\n');
fprintf('w_ci (Crossover): %.2f rad/s\n', wci);
fprintf('w_zi (Zero):      %.2f rad/s\n', wzi);
fprintf('tau_i:            %.4e s\n', tau_i);
fprintf('kc_i:             %.4f\n\n', kci);

fprintf('--- PARÂMETROS DA MALHA DE TENSÃO ---\n');
fprintf('w_cv (Crossover): %.2f rad/s\n', wcv);
fprintf('w_zv (Zero):      %.2f rad/s\n', omega_zv);
fprintf('tau_v:            %.4e s\n', tau_v);
fprintf('kc_v:             %.4f\n\n', kc_v);


% Malha Aberta e Fechada de Tensão (Cascata Completa)
FTLA_v = FTLA_ncv * Cv;
FTMF_v = feedback(FTLA_v, 1);

%% 4. Corrected Plotting Logic (The "IEEE Fix")
t = 0:1e-6:0.03; % 30ms simulation
blue = [0, 0.447, 0.741]; red = [0.85, 0.325, 0.098];

% --- Figure 1: Step Responses (Corrected Syntax) ---
figure('Color', 'w', 'Name', 'Step Responses');

% Current Step
[y_i, t_i] = step(FTMF_i, t);
subplot(2,1,1);
plot(t_i*1e3, y_i, 'LineWidth', 2, 'Color', blue); % LineWidth works here!
grid on; ylabel('Current (p.u.)'); title('Current Loop Step Response');

% Voltage Step
[y_v, t_v] = step(FTMF_v, t);
subplot(2,1,2);
plot(t_v*1e3, y_v, 'LineWidth', 2, 'Color', red);
grid on; ylabel('Voltage (p.u.)'); xlabel('Time (ms)');
title('Voltage Loop Step Response (Cascaded)');

% --- Figure 2: Bode Plots (Corrected Syntax) ---
figure('Color', 'w', 'Name', 'Bode Comparison');

% Extracting data for manual plotting to avoid 'LineWidth' errors in bode()
[magI, phaseI, wI] = bode(FTLA_i, {2*pi*10, 2*pi*20000});
[magV, phaseV, wV] = bode(FTLA_v, {2*pi*10, 2*pi*20000});

subplot(2,1,1);
semilogx(wI/(2*pi), 20*log10(squeeze(magI)), 'Color', blue, 'LineWidth', 2); hold on;
semilogx(wV/(2*pi), 20*log10(squeeze(magV)), 'Color', red, 'LineWidth', 2);
grid on; ylabel('Magnitude (dB)'); title('Bode Magnitude Comparison');
legend('Current Loop', 'Voltage Loop');

subplot(2,1,2);
semilogx(wI/(2*pi), squeeze(phaseI), 'Color', blue, 'LineWidth', 2); hold on;
semilogx(wV/(2*pi), squeeze(phaseV), 'Color', red, 'LineWidth', 2);
grid on; ylabel('Phase (deg)'); xlabel('Frequency (Hz)');
