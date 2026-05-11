%% controle corrente funcionando 
% Especificações da malha de tensão - planta lenta fase maior
wc_v = 5000 * 2 * pi;                % Frequência de cruzamento (rad/s) - 250 hz
MF_v = 65 * pi/180;         % Margem de fase desejada (90° em radianos)

% Função de transferência de malha aberta não compensada
FTLA_NC_v = FT_G_Vod * Kv ; 

% Calcular magnitude e fase na frequência de cruzamento
[mag_v, phase_v, ~] = bode(FTLA_NC_v, wc_v);
mag_v_abs = mag_v(1);
phase_v_deg = phase_v(1);
phase_v_rad = phase_v_deg * pi/180;
fprintf('Em ωc = %.2f rad/s (%.2f Hz):\n', wc_v, wc_v/(2*pi));
fprintf('  |FTLA_NC| = %.2f dB (%.4e linear)\n', 20*log10(mag_v_abs), mag_v_abs);
fprintf('  ∠FTLA_NC = %.2f° (%.4f rad)\n', phase_v_deg, phase_v_rad);

% Cálculo do zero do controlador (ωz)
tan_arg_v = tan(MF_v - pi/2 - phase_v_rad);
wz_v = wc_v /tan_arg_v;

tau_v = 1 / wz_v

kc_v = wc_v / (sqrt(wc_v^2 + wz_v^2) * mag_v_abs)
%%

% Especificações da malha de corrente
fci = 500;              
MF_i = 300 * pi/180;         % Margem de fase desejada (60° em radianos)
wci = 2*pi*fci  % Frequência de cruzamento (rad/s) - 2 kHz

% Função de transferência de malha aberta não compensada
% usando subs
FTLA_NC_i1 = FT_G_iLd * Ki * KPWM;

% usando ss
Gi_dc = Vin/(Dp^2*Rload);
FTLA_NC_i = FT_G_iL2d * Ki * KPWM;
[mag_i, phase_i, ~] = bode(FTLA_NC_i, wci);
mag_i_abs = mag_i(1);
phase_i_deg = phase_i(1);
phase_i_rad = phase_i_deg * pi/180;
fprintf('At ωc = %.2f rad/s (%.2f Hz):\n', wci, wci/(2*pi));
fprintf('  |FTLA_NC| = %.2f dB (%.4e linear)\n', 20*log10(mag_i_abs), mag_i_abs);
fprintf('  arg or ∠FTLA_NC = %.2f° (%.4f rad)\n', phase_i_deg, phase_i_rad);

% Cálculo do zero do controlador (ωz)
tan_arg_i = tan(MF_i - pi/2 - phase_i_rad);
wz_i = wci / tan_arg_i;
% Cálculo do ganho do controlador (kc)
kc_i = wci / (sqrt(wci^2 + wz_i^2) * mag_i_abs);
tau_i=1/wz_i;
fprintf('  tau = %.2f° (kc_i %.4f rad)\n', tau_i, kc_i);
