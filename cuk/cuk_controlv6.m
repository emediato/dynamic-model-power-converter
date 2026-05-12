
% Criar diretório para figuras se não existir
if ~exist('figures', 'dir'), mkdir('figures'); end

%% conversor

Vin   = 48;             % Tensão de entrada [V]
Vo    = 72;             % Tensão de saída desejada [V]
Po    = 300;            % Potência nominal [W]
fs    = 50e3;           % Frequência de chaveamento [Hz]
Rload = Vo^2 / Po;      % Resistência de carga nominal [Ohm]
D     = Vo / (Vin + Vo);% Razão cíclica teórica (Cuk: Vo/Vin = D/(1-D))
Kpwm  = 1.0;            % Ganho do modulador PWM

% Ganhos dos Sensores (Ajuste conforme seu hardware/simulação)
V_sensor_max = 3.3;
Ki = V_sensor_max / 20; % Sensor de corrente (5V para 20A)
Kv = V_sensor_max / Vo; % Sensor de tensão (5V para Vo nominal)

%% ESPECIFICAÇÕES DAS MALHAS (CRITÉRIO DE DESACOPLAMENTO)
% Malha de Corrente (Rápida)
fci  = 2000;            % 2.5 kHz (fs/20)
wci  = 2*pi*fci;
MF_i = 60;              % Margem de Fase desejada [deg]

% Malha de Tensão (Lenta)
fcv  = 250;             % 250 Hz (fci/10)
wcv  = 2*pi*fcv;
MF_v = 90;              % Margem de Fase maior para estabilidade (RHP zeros)

%% PROJETO DA MALHA INTERNA (CORRENTE iL2)
% Assume-se que GiL2d e GiL1d (TFs Corrente/Duty) já existem no workspace
GiL2d = sys_iL2;
Plant_iL2 = GiL2d * Ki * Kpwm;

[mag_i2, phase_i2] = bode(Plant_iL2, wci);
mag_i2 = squeeze(mag_i2); phase_i2 = squeeze(phase_i2);

% Síntese do Controlador PI_i (Malha iL2)
wz_i2 = wci / tan(deg2rad(MF_i) - pi/2 - deg2rad(phase_i2));
kc_i2 = wci / (mag_i2 * sqrt(wci^2 + wz_i2^2));

Ctrl_iL2 = kc_i2 * tf([1 wz_i2], [1 0]);
OL_iL2   = Plant_iL2 * Ctrl_iL2;
CL_iL2   = feedback(OL_iL2, 1); % Malha Fechada de Corrente

FTMF_i2 = Ctrl_iL2 * KPWM * FT_G_iL2d / (1 + Ki * Ctrl_iL2 * KPWM * FT_G_iL2d);


%% COMPARATIVO: MALHA DE CORRENTE iL1 (ENTRADA)
GiL1d = sys_iL1;
Plant_iL1 = GiL1d * Ki * Kpwm;

[mag_i1, phase_i1] = bode(Plant_iL1, wci);
mag_i1 = squeeze(mag_i1); phase_i1 = squeeze(phase_i1);

wz_i1 = wci / tan(deg2rad(MF_i) - pi/2 - deg2rad(phase_i1));
kc_i1 = wci / (mag_i1 * sqrt(wci^2 + wz_i1^2));

Ctrl_iL1 = kc_i1 * tf([1 wz_i1], [1 0]);
CL_iL1   = feedback(Plant_iL1 * Ctrl_iL1, 1);

%% PROJETO DA MALHA EXTERNA (TENSÃO Vo - CASCATA)
% A planta de tensão "vista" pelo controle é Vo / iL2_referência
% Gvi = (Vo/d) / (iL2/d)
Gvi_model = Gvd / GiL2d; 

% FTLA de Tensão considerando a malha de corrente fechada (Cascata)
Plant_Vo_Cascade = Gvd * Kv;

% Plant_Vo_Cascade = Gvi_model * FTMF_i2 * (Kv ) * Kpwm;

%Plant_Vo_Cascade = Gvi_model * FTMF_i2 * (Kv / Ki);


[mag_v, phase_v] = bode(Plant_Vo_Cascade, wcv);
mag_v = squeeze(mag_v); phase_v = squeeze(phase_v);

% Síntese do Controlador PI_v
wz_v = wcv / tan(deg2rad(MF_v) - pi/2 - deg2rad(phase_v));
kc_v = wcv / (mag_v * sqrt(wcv^2 + wz_v^2));

Ctrl_Vo = kc_v * tf([1 wz_v], [1 0]);
OL_Vo   = Plant_Vo_Cascade * Ctrl_Vo;
CL_Vo   = feedback(OL_Vo, 1);

%% GERAÇÃO DE GRÁFICOS E FIGURAS COMPARATIVAS

% --- Figura 1: Comparativo de Resposta de Corrente (iL1 vs iL2) ---
fig1 = figure('Name', 'Comparativo Malhas de Corrente', 'Color', 'w');
subplot(2,1,1);
[y1, t1] = step(CL_iL1, 0.005); 
[y2, t2] = step(FTMF_i2, 0.005);
plot(t1*1e3, y1, 'LineWidth', 1.5); hold on;
plot(t2*1e3, y2, 'LineWidth', 1.5);
grid on; ylabel('Amp. Normalizada'); title('Resposta ao Degrau: Malhas de Corrente');
legend('iL1 (Input)', 'iL2 (Output)');

subplot(2,1,2);
bode(OL_iL2); hold on; 
% bode(OL_iL1);
grid on; title('Bode: Malha Aberta Compensada (Corrente)');
print(fig1, 'figures/comp_current_loops.png', '-dpng', '-r300');

% --- Figura 2: Malha de Tensão em Cascata ---
fig2 = figure('Name', 'Resultado Malha de Tensão', 'Color', 'w');
subplot(2,1,1);
step(CL_Vo, 0.05); grid on;
title('Resposta ao Degrau: Malha de Tensão (Cascata Completa)');
ylabel('Tensão (p.u.)');

subplot(2,1,2);
margin(OL_Vo); grid on;
print(fig2, 'figures/voltage_cascade_response.png', '-dpng', '-r300');

% --- Figura 3: Comparativo State Space vs Transfer Function ---
% (Apenas se as variáveis SS existirem)
if exist('FTLA_NC_v_ss', 'var')
    fig3 = figure('Name', 'Comparativo SS vs TF', 'Color', 'w');
    bode(Plant_Vo_Cascade, 'b'); hold on;
    bode(FTLA_NC_v_ss, 'r--');
    grid on; title('Comparativo: Modelo SS vs Modelo TF');
    legend('Transfer Function', 'State Space');
    print(fig3, 'figures/comp_ss_vs_tf.png', '-dpng', '-r300');
end

%%  EXIBIÇÃO DE PARÂMETROS FINAIS
fprintf('\n==========================================\n');
fprintf('   RELATÓRIO DE SINTESE DE CONTROLE\n');
fprintf('==========================================\n');
 1/wz_i2
fprintf('PI CORRENTE (iL2): Kc = %.4e | tau = %.2f | wz = %.2f (fz = %.1f Hz)\n', kc_i2, 1/wz_i2, wz_i2, wz_i2/(2*pi));
fprintf('PI TENSÃO (Vo):   Kc = %.4e | tau = %.2f | wz = %.2f (fz = %.1f Hz)\n', kc_v, 1/wz_v, wz_v, wz_v/(2*pi));
1/wz_v


fprintf('------------------------------------------\n');
[~, pm_v_final] = margin(OL_Vo);
fprintf('Margem de Fase Final (Tensão): %.2f graus\n', pm_v_final);
fprintf('Frequência de Cruzamento: %.2f Hz\n', wcv/(2*pi));
fprintf('==========================================\n');