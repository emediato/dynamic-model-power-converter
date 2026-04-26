%% Script: Diagrama de Bode para Inversor Monofásico
% Este script calcula e plota a magnitude e fase da FTLA (Função de 
% Transferência de Malha Aberta) não compensada e compensada
% com base nos parâmetros fornecidos.

clear; clc; close all;

%% 1. Parâmetros do sistema (conforme fornecido)
Vi = 400;          % Tensão de entrada (V)
Vo = 311;          % Tensão de saída (V)
Po = 1250;         % Potência de saída (W)
Ro = Vo^2 / Po;    % Resistência de carga (Ohm) ~ 77.38 Ohm
Io = Vo / Ro;      % Corrente de saída (A) ~ 4.02 A
D = Vo / Vi;       % Ciclo de trabalho ~ 0.778
fs = 20e3;         % Frequência de comutação (Hz)
Lo = 8e-3;         % Indutância do filtro (H) = 8 mH
Co = 1e-6;         % Capacitância do filtro (F) = 1 uF
j = sqrt(-1);      % Unidade imaginária
kv = 3.3 / 311;    % Ganho do sensor de tensão (V/V)
kpwm = 1 / 1;      % Ganho do modulador PWM (1/V)

%% 2. Parâmetros do compensador (conforme fornecido)
% Frequência de cruzamento desejada (1/60 da frequência de comutação)
wc = 2 * pi * fs / 60;  % rad/s
% Margem de fase desejada (90 graus em radianos)
Mphi = 90 * pi / 180;   % rad

%% 3. Definição da função de transferência não compensada: G(s)
% G(s) = Vi * (1/(Lo*Co)) * [s^2 + (1/(Ro*Co)) * s + 1/(Lo*Co)]^-1
% Ou forma mais clara: G(s) = Vi * (1/(Lo*Co)) / (s^2 + (1/(Ro*Co))*s + 1/(Lo*Co))

% Coeficientes do denominador
a1 = 1/(Ro*Co);
a0 = 1/(Lo*Co);

% Numerador: Vi * a0
num_G = Vi * a0;
den_G = [1, a1, a0];

% Criação da função de transferência G(s)
s = tf('s');
G = Vi * (1/(Lo*Co)) / (s^2 + (1/(Ro*Co))*s + 1/(Lo*Co));

%% 4. FTLA não compensada: FTLA_nc(s) = kv * kpwm * G(s)
kv_kpwm = kv * kpwm;
FTLA_nc = kv_kpwm * G;

%% 5. Cálculo do zero do compensador (omega_z) e ganho kc
% Primeiro, obtemos a fase da FTLA não compensada na frequência wc
[mag_nc, phase_nc] = bode(FTLA_nc, wc);
mag_nc = squeeze(mag_nc);    % magnitude em valor linear
phase_nc = squeeze(phase_nc); % fase em graus

% Converter fase para radianos
phase_nc_rad = phase_nc * pi/180;

% Cálculo de omega_z conforme a equação fornecida
% omega_z = wc / tan(Mphi - pi/2 - arg(FTLA_nc(wc)))
den_tan = tan(Mphi - pi/2 - phase_nc_rad);
% Evitar divisão por zero
if abs(den_tan) < 1e-6
    den_tan = 1e-6;
end
omega_z = wc / den_tan;   % rad/s
tau = 1/omega_z;          % constante de tempo (s)
tau_us = tau * 1e6;       % em microssegundos

% Cálculo de kc
% kc = wc / ( |FTLA_nc(wc)| * sqrt(wc^2 + omega_z^2) )
kc = wc / (mag_nc * sqrt(wc^2 + omega_z^2));

% Exibir resultados dos parâmetros do compensador
fprintf('=== Parâmetros do Compensador ===\n');
fprintf('Frequência de cruzamento (wc) = %.2f rad/s (%.2f Hz)\n', wc, wc/(2*pi));
fprintf('omega_z = %.2f rad/s\n', omega_z);
fprintf('tau = %.2e s = %.2f µs\n', tau, tau_us);
fprintf('kc = %.4f\n\n', kc);

%% 6. Definição do compensador C(s)
% C(s) = kc * (s + omega_z) / s
C = kc * (s + omega_z) / s;

%% 7. FTLA compensada: FTLA_c(s) = C(s) * FTLA_nc(s)
FTLA_c = C * FTLA_nc;

%% 8. Diagramas de Bode (Magnitude e Fase)
% Faixa de frequências para o gráfico (1 Hz a 100 kHz)
freq_min = 1;      % Hz
freq_max = 100e3;  % Hz
frequencies = logspace(log10(freq_min), log10(freq_max), 500);
w = 2*pi*frequencies; % rad/s

% Avaliar resposta em frequência
[mag_nc_plot, phase_nc_plot] = bode(FTLA_nc, w);
[mag_c_plot, phase_c_plot] = bode(FTLA_c, w);

% Remover dimensões extras
mag_nc_plot = squeeze(mag_nc_plot);
phase_nc_plot = squeeze(phase_nc_plot);
mag_c_plot = squeeze(mag_c_plot);
phase_c_plot = squeeze(phase_c_plot);

% Converter magnitude para dB
mag_nc_dB = 20*log10(mag_nc_plot);
mag_c_dB = 20*log10(mag_c_plot);

%% 9. Plotagem dos diagramas de Bode
figure('Name', 'Diagrama de Bode do Inversor Monofásico', 'Position', [100 100 900 600]);

% Subplot 1: Magnitude
subplot(2,1,1);
semilogx(frequencies, mag_nc_dB, 'b-', 'LineWidth', 1.5); hold on;
semilogx(frequencies, mag_c_dB, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
title('Diagrama de Bode - Magnitude');
legend('FTLA não compensada', 'FTLA compensada', 'Location', 'southwest');
xlim([freq_min, freq_max]);

% Subplot 2: Fase
subplot(2,1,2);
semilogx(frequencies, phase_nc_plot, 'b-', 'LineWidth', 1.5); hold on;
semilogx(frequencies, phase_c_plot, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Diagrama de Bode - Fase');
legend('FTLA não compensada', 'FTLA compensada', 'Location', 'southwest');
xlim([freq_min, freq_max]);

% Adicionar anotação na frequência de cruzamento
subplot(2,1,1);
line([wc/(2*pi) wc/(2*pi)], [-100 100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text(wc/(2*pi), -80, sprintf('f_c = %.1f Hz', wc/(2*pi)), 'HorizontalAlignment', 'center');

subplot(2,1,2);
line([wc/(2*pi) wc/(2*pi)], [-200 100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

%% 10. Verificação da margem de fase na frequência de cruzamento
[mag_c_c, phase_c_c] = bode(FTLA_c, wc);
mag_c_c = squeeze(mag_c_c);
phase_c_c = squeeze(phase_c_c);
margem_fase = 180 + phase_c_c; % margem de fase em graus

fprintf('=== Verificação do Projeto ===\n');
fprintf('Na frequência de cruzamento f_c = %.2f Hz:\n', wc/(2*pi));
fprintf('Magnitude de FTLA_c = %.2f dB\n', 20*log10(mag_c_c));
fprintf('Fase de FTLA_c = %.2f graus\n', phase_c_c);
fprintf('Margem de fase obtida = %.2f graus (desejada: %.2f graus)\n', margem_fase, Mphi*180/pi);

if abs(margem_fase - Mphi*180/pi) < 10
    fprintf('>>> Projeto atende à margem de fase desejada.\n');
else
    fprintf('>>> Ajuste fino do compensador pode ser necessário.\n');
end

%% 11. Funções de transferência em forma textual
fprintf('\n=== Função de Transferência ===\n');
fprintf('G(s) = \n');
pretty(G)
fprintf('C(s) = \n');
pretty(C)
fprintf('FTLA compensada = \n');
pretty(FTLA_c)
