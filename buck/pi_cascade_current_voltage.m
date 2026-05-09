%% ========================================================================
%  Cascade Control of Buck Converter
%  Inner: Current Loop (already designed)
%  Outer: Voltage Loop (designed here)
%  IEEE-style plots for publication
% ========================================================================

clear; close all; clc;

%% 1. System Parameters (from previous images)
Vin = 100;      % V
Ro = 6;         % Ohm
Lo = 2e-3;      % H
Co = 10e-6;     % F
fs = 40e3;      % Hz
ki = 0.33;      % V/A (current sensor gain)
kv = 0.055;     % voltage sensor gain (3.3V/60V)
kpwm = 1;       % modulator gain

%% 2. Current Loop Design (from previous data)
s = tf('s');
% Plant current: Gid = Vin / (s*Lo)
Gid = Vin / (s * Lo);

% Current compensator (lead) parameters from previous design
wci = 2*pi*2000;    % crossover frequency 2kHz
wzi = 3152;         % zero frequency (rad/s)
kci = 1.758;        % compensator gain
Ci = kci * (s + wzi) / s;

% Current open-loop transfer function
FTLA_nci = Gid * ki * kpwm;
FTLA_ci = FTLA_nci * Ci;

% Current closed-loop transfer function
FTMF_i = feedback(FTLA_ci, 1);

% Verify current loop performance
[Gm_i, Pm_i, ~, ~] = margin(FTLA_ci);
fprintf('=== Current Loop Performance ===\n');
fprintf('Crossover frequency: %.1f Hz\n', wci/(2*pi));
fprintf('Phase margin: %.2f deg\n', Pm_i);
fprintf('Gain margin: %.2f dB\n', 20*log10(Gm_i));
fprintf('================================\n\n');

%% 3. Voltage Plant (Gvi)
% Gvi(s) = Ro / (s*Ro*Co + 1)
Gvi = Ro / (s * Ro * Co + 1);

%% 4. Voltage Open-Loop without compensator
% FTLA_ncv = Gvi * kv * FTMF_i
FTLA_ncv = Gvi * kv * FTMF_i;

%% 5. Voltage Loop Specifications
fcv = 200;          % Hz (crossover frequency)
wcv = 2*pi*fcv;     % rad/s
Mphiv = 90 * pi/180; % desired phase margin (90°)

% Phase of FTLA_ncv at wcv
[mag_ncv, phase_ncv] = bode(FTLA_ncv, wcv);
phase_ncv_deg = phase_ncv;

% Required phase lead
phi_lead_v = rad2deg(Mphiv) - 90 - phase_ncv_deg;
phi_lead_v_rad = deg2rad(phi_lead_v);

% Zero frequency of voltage compensator
wzv = wcv / tan(phi_lead_v_rad + pi/2);
tau_v = 1/wzv;      % time constant (s)

% Gain at crossover frequency
kcv = wcv / (mag_ncv * sqrt(wzv^2 + wcv^2));

%% 6. Voltage Compensator (lead)
Cv = kcv * (s + wzv) / s;

%% 7. Voltage Open-Loop Compensated
FTLA_cv = FTLA_ncv * Cv;

%% 8. Voltage Closed-Loop Transfer Function
FTMF_v = feedback(FTLA_cv, 1);

%% 9. Performance Analysis
[Gm_v, Pm_v, ~, ~] = margin(FTLA_cv);

fprintf('=== Voltage Loop Performance ===\n');
fprintf('Crossover frequency: %.1f Hz\n', fcv);
fprintf('Phase margin: %.2f deg\n', Pm_v);
fprintf('Gain margin: %.2f dB\n', 20*log10(Gm_v));
fprintf('Zero frequency wzv: %.2f rad/s (tau = %.2f us)\n', wzv, tau_v*1e6);
fprintf('Compensator gain kcv: %.3f\n', kcv);
fprintf('================================\n\n');

%% 10. Time-Domain Simulation
t = 0:1e-5:5e-3;    % 5ms simulation
% Step response (output voltage)
[Vout, t_out] = step(FTMF_v, t);
% Step response (inductor current)
[I_L, t_i] = step(FTMF_i, t);

%% 11. IEEE-Style Plots
% Color scheme for publication
blue = [0, 0.4470, 0.7410];
red = [0.8500, 0.3250, 0.0980];
green = [0.4660, 0.6740, 0.1880];

% -------------------------------------------------
% Figure 1: Bode Plots - Current and Voltage Loops
% -------------------------------------------------
figure('Color', 'white', 'Position', [100, 100, 900, 700]);

% Current loop Bode
subplot(2,1,1);
bode(FTLA_ci, {2*pi*100, 2*pi*20000}, 'LineWidth', 1.5, 'Color', blue);
grid on;
title('(a) Current Loop - Open-Loop Transfer Function', 'FontSize', 11);
legend('FTLA_{ci}(s)', 'Location', 'southwest');

% Voltage loop Bode
subplot(2,1,2);
bode(FTLA_cv, {2*pi*10, 2*pi*2000}, 'LineWidth', 1.5, 'Color', red);
grid on;
title('(b) Voltage Loop - Open-Loop Transfer Function', 'FontSize', 11);
legend('FTLA_{cv}(s)', 'Location', 'southwest');

% -------------------------------------------------
% Figure 2: Step Responses (Current and Voltage)
% -------------------------------------------------
figure('Color', 'white', 'Position', [100, 100, 900, 700]);

subplot(2,1,1);
plot(t_i*1e3, I_L, 'LineWidth', 1.5, 'Color', blue);
grid on;
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Inductor Current (A)', 'FontSize', 11);
title('(a) Current Loop Step Response', 'FontSize', 11);
axis([0 2 0 1.2]);
hold on;
plot([0 2], [1 1], 'k--', 'LineWidth', 0.8);
legend('I_L(t)', 'Reference', 'Location', 'southeast');

subplot(2,1,2);
plot(t_out*1e3, Vout, 'LineWidth', 1.5, 'Color', red);
grid on;
xlabel('Time (ms)', 'FontSize', 11);
ylabel('Output Voltage (V)', 'FontSize', 11);
title('(b) Voltage Loop Step Response', 'FontSize', 11);
axis([0 5 0 1.2]);
hold on;
plot([0 5], [1 1], 'k--', 'LineWidth', 0.8);
legend('V_o(t)', 'Reference', 'Location', 'southeast');

% -------------------------------------------------
% Figure 3: Combined Bode (Comparison)
% -------------------------------------------------
figure('Color', 'white', 'Position', [100, 100, 900, 500]);
bode(FTLA_ci, FTLA_cv, {2*pi*10, 2*pi*20000});
grid on;
title('Open-Loop Bode Comparison: Current vs. Voltage Loops', 'FontSize', 12);
legend('Current Loop (FTLA_{ci})', 'Voltage Loop (FTLA_{cv})', 'Location', 'southwest');

% -------------------------------------------------
% Figure 4: Cascaded Response (Output Voltage)
% -------------------------------------------------
figure('Color', 'white', 'Position', [100, 100, 600, 400]);
plot(t_out*1e3, Vout, 'LineWidth', 2, 'Color', green);
grid on;
xlabel('Time (ms)', 'FontSize', 12);
ylabel('Output Voltage (p.u.)', 'FontSize', 12);
title('Cascaded Voltage Control - Step Response', 'FontSize', 12);
axis([0 5 0 1.2]);
hold on;
plot([0 5], [1 1], 'k--', 'LineWidth', 1);
text(3.5, 0.85, ['Settling time: ', num2str(t_out(find(Vout>0.98,1,'first'))*1000, '%.1f'), ' ms'], 'FontSize', 10);
text(3.5, 0.75, ['Overshoot: ', num2str((max(Vout)-1)*100, '%.1f'), ' %'], 'FontSize', 10);

%% 11. Additional Metrics for Publication
fprintf('\n=== Performance Metrics for Publication ===\n');
% Settling time (2%)
idx_settle = find(abs(Vout-1) < 0.02, 1, 'first');
settling_time = t_out(idx_settle)*1000;
overshoot = (max(Vout)-1)*100;
rise_time = t_out(find(Vout>0.9,1,'first'))*1000;

fprintf('Settling time (2%%): %.2f ms\n', settling_time);
fprintf('Overshoot: %.2f %%\n', overshoot);
fprintf('Rise time (10-90%%): %.2f ms\n', rise_time);
fprintf('Current loop bandwidth: %.1f Hz\n', wci/(2*pi));
fprintf('Voltage loop bandwidth: %.1f Hz\n', fcv);
fprintf('============================================\n');

%% 12. Export Data for Publication (optional)
% Uncomment to save figures as EPS or PNG
% saveas(gcf, 'voltage_step_response.eps', 'epsc');
% saveas(gcf, 'bode_comparison.eps', 'epsc');
