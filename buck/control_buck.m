%% Controle do Conversor Buck com PI por Resposta em Frequência
% Parâmetros fornecidos
Vi = 48;         % Tensão de entrada (V)
Vo = 24;         % Tensão de saída (V)
Po = 96;         % Potência de saída (W)
Ro = Vo^2 / Po;  % Resistência de carga (Ohm) = 6
D = Vo / Vi;     % Ciclo de trabalho = 0.5
fs = 20e3;       % Frequência de comutação (Hz)
Lo = 3e-3;       % Indutância (H)
Co = 10e-6;      % Capacitância (F)

% Ganhos do sensor e PWM
kv = 4 / 24;          % Ganho do sensor de tensão (4V -> 24V)
kpwm = 1 / 5;         % Ganho do PWM (pico de 5V)

% Unidade imaginária
j = sqrt(-1);

% Frequência de cruzamento desejada (1/20 da frequência de comutação)
wc = 2 * pi * fs / 20;  % rad/s
fprintf('Frequência de cruzamento desejada: %.2f rad/s (%.2f Hz)\n', wc, wc/(2*pi));

% Margem de fase desejada (60 graus)
Mphi = 60 * pi / 180;

%% 1. DEFINIÇÃO DA PLANTA G(s)
s = tf('s');
G = Vi / (Lo*Co*s^2 + (Lo/Ro)*s + 1);

% Verificação da função
fprintf('\n--- Função de Transferência da Planta G(s) ---\n');
G

%% 2. FTLA NÃO COMPENSADA (FTLA_nc)
FTLA_nc = kv * kpwm * G;

fprintf('\n--- FTLA não compensada ---\n');
FTLA_nc

%% 3. PROJETO DO COMPENSADOR PI
% O compensador PI tem a forma: C(s) = kc * (s + wz) / s
% Onde wz é o zero do compensador e kc é o ganho.

% Primeiro, obter a fase da FTLA_nc na frequência wc
[mag_nc, phase_nc] = bode(FTLA_nc, wc);
mag_nc = squeeze(mag_nc);
phase_nc_deg = squeeze(phase_nc);
phase_nc_rad = phase_nc_deg * pi / 180;

% Cálculo de wz (zero do compensador)
den_tan = tan(Mphi - pi/2 - phase_nc_rad);
if abs(den_tan) < 1e-6
    den_tan = 1e-6; % Evitar divisão por zero
end
wz = wc / den_tan;
tau = 1 / wz;  % Constante de tempo (s)
tau_us = tau * 1e6;  % em microssegundos

% Cálculo de kc (ganho do compensador)
kc = wc / (mag_nc * sqrt(wc^2 + wz^2));

fprintf('\n--- Parâmetros do Compensador PI ---\n');
fprintf('wz (zero do compensador) = %.2f rad/s\n', wz);
fprintf('Constante de tempo tau = %.2e s = %.2f µs\n', tau, tau_us);
fprintf('kc (ganho do compensador) = %.4f\n', kc);

% Compensador PI
C = kc * (s + wz) / s;
fprintf('\n--- Compensador C(s) ---\n');
C

%% 4. FTLA COMPENSADA
FTLA_c = C * FTLA_nc;

%% 5. DIAGRAMA DE BODE
freq_min = 1;      % Hz
freq_max = 100e3;  % Hz
frequencies = logspace(log10(freq_min), log10(freq_max), 500);
w = 2*pi*frequencies;

% Avaliar resposta
[mag_nc_plot, phase_nc_plot] = bode(FTLA_nc, w);
[mag_c_plot, phase_c_plot] = bode(FTLA_c, w);

mag_nc_plot = squeeze(mag_nc_plot);
phase_nc_plot = squeeze(phase_nc_plot);
mag_c_plot = squeeze(mag_c_plot);
phase_c_plot = squeeze(phase_c_plot);

mag_nc_dB = 20*log10(mag_nc_plot);
mag_c_dB = 20*log10(mag_c_plot);

% Plotagem
figure('Name', 'Diagrama de Bode do Conversor Buck', 'Position', [100 100 900 600]);

subplot(2,1,1);
semilogx(frequencies, mag_nc_dB, 'b-', 'LineWidth', 1.5); hold on;
semilogx(frequencies, mag_c_dB, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
title('Diagrama de Bode - Magnitude');
legend('FTLA não compensada', 'FTLA compensada', 'Location', 'southwest');
xlim([freq_min, freq_max]);
line([wc/(2*pi) wc/(2*pi)], [-100 100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
text(wc/(2*pi), -80, sprintf('f_c = %.1f Hz', wc/(2*pi)), 'HorizontalAlignment', 'center');

subplot(2,1,2);
semilogx(frequencies, phase_nc_plot, 'b-', 'LineWidth', 1.5); hold on;
semilogx(frequencies, phase_c_plot, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('Frequência (Hz)');
ylabel('Fase (graus)');
title('Diagrama de Bode - Fase');
legend('FTLA não compensada', 'FTLA compensada', 'Location', 'southwest');
xlim([freq_min, freq_max]);
line([wc/(2*pi) wc/(2*pi)], [-200 100], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);

%% 6. MARGENS DE FASE E GANHO
% Margens do sistema compensado
[Gm, Pm, Wcg, Wcp] = margin(FTLA_c);

fprintf('\n--- Margens do Sistema Compensado ---\n');
fprintf('Margem de ganho (Gm) = %.2f dB\n', 20*log10(Gm));
fprintf('Margem de fase (Pm) = %.2f graus\n', Pm);
fprintf('Frequência de cruzamento de ganho (Wcp) = %.2f rad/s (%.2f Hz)\n', Wcp, Wcp/(2*pi));

% Verificar se atende especificação
if Pm >= Mphi*180/pi - 5
    fprintf('>>> Projeto atende à margem de fase desejada (%.1f graus).\n', Mphi*180/pi);
else
    fprintf('>>> Margem de fase obtida (%.1f) < desejada (%.1f). Revisar projeto.\n', Pm, Mphi*180/pi);
end

% Plotar margens na figura de Bode
figure(margin(FTLA_c));
title('Margens de Fase e Ganho do Sistema Compensado');
