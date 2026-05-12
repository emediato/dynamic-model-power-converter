cuk_controlv5
   # RELATÓRIO DE SINTESE DE CONTROLE
### tensao
wc_v = 200 * 2 * pi;                % Frequência de cruzamento (rad/s) - 250 hz
MF_v = 30 * pi/180;         % Margem de fase desejada (90° em radianos)
Em ωc = 1256.64 rad/s (200.00 Hz):
  |FTLA_NC| = 22.80 dB (1.3798e+01 linear)
  ∠FTLA_NC = 354.00° (6.1784 rad)
kc_v =    0.0586
tau_v =   -0.0011
### corrente % Especificações da malha de corrente
fci = 2000;              
MF_i = 90 * pi/180;         % Margem de fase desejada (60° em radianos)
wci = 2*pi*fci  % Frequência de cruzamento (rad/s) - 2 kHz

At ωc = 12566.37 rad/s (2000.00 Hz):
  |FTLA_NC| = 59.33 dB (9.2573e+02 linear)
  arg or ∠FTLA_NC = -5.72° (-0.0998 rad)
kc_i =  1.0765e-04
tau_i =   7.9701e-06


cuk_controlv6
   # RELATÓRIO DE SINTESE DE CONTROLE
fci  = 2000;            % 2.5 kHz (fs/20)
MF_i = 90;              % Margem de Fase desejada [deg]

% Malha de Tensão (Lenta)
fcv  = 200;             % 250 Hz (fci/10)
MF_v = 55;              % Margem de Fase maior para estabilidade (RHP zeros)
PI CORRENTE (iL2): Kc = 4.8686e-01 | tau = 0.0030 | wz = 330.33 (fz = 52.6 Hz)
PI TENSÃO (Vo):   Kc = 9.1738e-03 | tau = 1.0155e-04 | wz = 9847.77 (fz = 1567.3 Hz)

   
