 % --- Configuração de Vetores e Pré-alocação ---
fci_vector  = [160, 5e3, 4e3, 1e3, 2e3, 200, 2.5e3];
mf_i_vector = [290, 65, 60, 60, 70, 60, 60];
fcv_vector  = [16, 500, 400, 100, 2e3, 20, 250]; 
mf_v_vector = [120, 300, 90, 70, 90, 90, 90];


n = length(fci_vector);
kci_vector = zeros(1, n);
tau_ci_vector = zeros(1, n);
kcv_vector = zeros(1, n); 
tau_v_vector = zeros(1, n);

% Celular para armazenar as TFs (Malha Aberta e Fechada)
FTLA_ci_cell = cell(1, n); CL_ci_cell = cell(1, n);
FTLA_cv_cell = cell(1, n); CL_cv_cell = cell(1, n);

for i = 1:n
    % --- 1. Projeto do Controlador de Corrente (Interno) ---
    wci = 2 * pi * fci_vector(i);
    [mag_i, phase_i] = bode(FTLA_NC_i, wci);
    phase_i_rad = phase_i(1) * pi/180;
    
    % Sintonia PI de Corrente
    wz_i = wci / tan((mf_i_vector(i)*pi/180) - pi/2 - phase_i_rad);
    if wz_i <= 0, wz_i = wci/10; end
    kc_i = wci / (mag_i(1) * sqrt(wci^2 + wz_i^2));
    
    % Armazenamento de dados
    kci_vector(i) = kc_i; tau_ci_vector(i) = 1/wz_i;
    Ci = tf(kc_i * [1, wz_i], [1, 0]);
    FTLA_ci_cell{i} = FTLA_NC_i * Ci;
    CL_ci_cell{i} = feedback(FTLA_ci_cell{i}, 1); % Malha Fechada
    
    FTMF_i = Ci * Kpwm * FT_G_iLd / ( 1 + Ki * (Ci * Kpwm * FT_G_iLd)) ;

    % --- 2. Projeto do Controlador de Tensão (Externo) ---
    wcv = 2 * pi * fcv_vector(i);
    FTLA_NC_v = FT_G_Vi * Kv * FTMF_i;   % Malha aberta não compensada
    [mag_v, phase_v] = bode(FTLA_NC_v, wcv);
    phase_v_rad = phase_v(1) * pi/180;
    
    % Sintonia PI de Tensão
    wz_v = wcv / tan((mf_v_vector(i)*pi/180) - pi/2 - phase_v_rad);
    if wz_v <= 0, wz_v = wcv/10; end
    kc_v = wcv / (mag_v(1) * sqrt(wcv^2 + wz_v^2));
    
    % Armazenamento de dados
    kcv_vector(i) = kc_v; tau_v_vector(i) = 1/wz_v;
    Cv = tf(kc_v * [1, wz_v], [1, 0]);
    FTLA_cv_cell{i} = FTLA_NC_v * Cv;
    CL_cv_cell{i} = feedback(FTLA_cv_cell{i}, 1); % Malha Fechada
end

% --- PLOTS COMPARATIVOS E ORGANIZAÇÃO DE DADOS ---

% Criando strings de legenda com Kc e Tau formatados
leg_i = arrayfun(@(x) sprintf('Fci=%dHz | Kc=%.2e | Tau=%.2e', ...
    fci_vector(x), kci_vector(x), tau_ci_vector(x)), 1:n, 'UniformOutput', false);

leg_v = arrayfun(@(x) sprintf('Fcv=%dHz | Kc=%.2e | Tau=%.2e', ...
    fcv_vector(x), kcv_vector(x), tau_v_vector(x)), 1:n, 'UniformOutput', false);


% FIGURE 1: CURRENT LOOP ANALYSIS (Bode and Step Response)
 
figure('Name', 'Current Loop Comparison', 'Units', 'normalized', 'Position', [0.1, 0.1, 0.4, 0.8]);

% Subplot 1: Bode Diagram (Frequency Response)
subplot(2,1,1); hold on;
for i = 1:n
    bode(FTLA_ci_cell{i});
end
grid on;
title('Frequency Response (Bode) - Current Loop');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB) / Phase (deg)');
hold off;

% Subplot 2: Step Response (Time Response)
subplot(2,1,2); hold on;
for i = 1:n
    step(CL_ci_cell{i});
end
grid on;
title('Step Response - Current Loop');
xlabel('Time (seconds)');
ylabel('Amplitude');
legend(leg_i, 'FontSize', 8);
hold off;


% FIGURE 2: VOLTAGE LOOP ANALYSIS (Bode and Step Response)

% This figure compares the frequency and time responses of the voltage loop
% for different controller configurations (n variations). 
figure('Name', 'Voltage Loop Comparison', 'Units', 'normalized', 'Position', [0.5, 0.1, 0.4, 0.8]);

% Subplot 1: Bode Diagram (Frequency Response)
subplot(2,1,1); hold on;
for i = 1:n
    bode(FTLA_cv_cell{i});
end
grid on;
title('Frequency Response (Bode) - Voltage Loop');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB) / Phase (deg)');
hold off;

% Subplot 2: Step Response (Time Response)
subplot(2,1,2); hold on;
for i = 1:n
    step(CL_cv_cell{i});
end
grid on;
title('Step Response - Voltage Loop');
xlabel('Time (seconds)');
ylabel('Amplitude');
legend(leg_v, 'FontSize', 8);
hold off;
