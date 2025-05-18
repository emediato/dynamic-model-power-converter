function [L1, L2] = calcular_indutores_cuk(Vi, Vo, Po, D, fs, delta_IL1_percent, delta_IL2_percent)
    % Função para calcular as indutâncias L1 e L2 do conversor Cuk
    %
    % Parâmetros de entrada:
    %   Vi: Tensão de entrada (V)
    %   Vo: Tensão de saída (V)
    %   Po: Potência de saída (W)
    %   D: Razão cíclica (0 < D < 1)
    %   fs: Frequência de comutação (Hz)
    %   delta_IL1_percent: Ondulação percentual da corrente em L1 (%)
    %   delta_IL2_percent: Ondulação percentual da corrente em L2 (%)
    %
    % Saídas:
    %   L1: Indutância L1 (H)
    %   L2: Indutância L2 (H)
    
    % Cálculo das correntes médias
    I_L1 = Po / Vi;      % Corrente média em L1
    I_L2 = Po / Vo;      % Corrente média em L2
    
    % Cálculo das ondulações de corrente
    delta_IL1 = I_L1 * (delta_IL1_percent/100);  % Ondulação absoluta em L1
    delta_IL2 = I_L2 * (delta_IL2_percent/100);  % Ondulação absoluta em L2
    
    % Cálculo das indutâncias
    L1 = (Vi * D) / (delta_IL1 * fs);  % Indutância L1
    L2 = (Vi * D) / (delta_IL2 * fs);  % Indutância L2
    
    % Exibição dos resultados
    fprintf('--- Resultados do Cálculo para Conversor Cuk ---\n');
    fprintf('Corrente média em L1 (I_L1): %.4f A\n', I_L1);
    fprintf('Corrente média em L2 (I_L2): %.4f A\n', I_L2);
    fprintf('Ondulação de corrente em L1 (ΔIL1): %.4f A\n', delta_IL1);
    fprintf('Ondulação de corrente em L2 (ΔIL2): %.4f A\n', delta_IL2);
    fprintf('Indutância L1 calculada: %.6f H (%.2f μH)\n', L1, L1*1e6);
    fprintf('Indutância L2 calculada: %.6f H (%.2f μH)\n', L2, L2*1e6);
end

[L1, L2] = calcular_indutores_cuk(12, 5, 50, 0.4, 100e3, 20, 15);


function [C1, C2] = calcular_capacitores_cuk(Vi, Vo, Io, D, fs, delta_Vc1_percent, delta_Vo_percent)
    % Função para calcular as capacitâncias C1 e C2 do conversor Cuk
    %
    % Parâmetros de entrada:
    %   Vi: Tensão de entrada (V)
    %   Vo: Tensão de saída (V)
    %   Io: Corrente de saída (A)
    %   D: Razão cíclica (0 < D < 1)
    %   fs: Frequência de comutação (Hz)
    %   delta_Vc1_percent: Ondulação percentual da tensão em C1 (%)
    %   delta_Vo_percent: Ondulação percentual da tensão de saída (%)
    %
    % Saídas:
    %   C1: Capacitância C1 (F)
    %   C2: Capacitância C2 (F)
    
    % Cálculo das ondulações de tensão
    delta_Vc1 = Vi * (delta_Vc1_percent/100);  % Ondulação absoluta em C1
    delta_Vo = Vo * (delta_Vo_percent/100);    % Ondulação absoluta na saída
    
    % Cálculo das capacitâncias
    C1 = (Io * D) / (delta_Vc1 * fs);  % Capacitância C1
    C2 = (Io * D) / (delta_Vo * fs);   % Capacitância C2
    
    % Exibição dos resultados
    fprintf('--- Resultados do Cálculo para Conversor Cuk ---\n');
    fprintf('Ondulação de tensão em C1 (ΔVc1): %.4f V\n', delta_Vc1);
    fprintf('Ondulação de tensão de saída (ΔVo): %.4f V\n', delta_Vo);
    fprintf('Capacitância C1 calculada: %.6f F (%.2f μF)\n', C1, C1*1e6);
    fprintf('Capacitância C2 calculada: %.6f F (%.2f μF)\n', C2, C2*1e6);
end

Tensão de entrada (Vi): 24V
Tensão de saída (Vo): 12V
Corrente de saída (Io): 2A
Razão cíclica (D): 0.5
Frequência de comutação (fs): 100kHz
Ondulação em C1: 5%
Ondulação na saída: 2%