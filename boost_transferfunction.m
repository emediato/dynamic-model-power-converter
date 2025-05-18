Vi = 12
Vo = 23
D = 0.5

L = 144*10e-6
C = 17.36*10E-6
Ro =11.52

j = sqrt(-1);
% x1 = linspace(1, 2, 5)
% varia de 1 em 1 até 10 a 5

% omega = [1..2] % controle em baixa frequencia

%numerator = Vi/C*L - Vi*s/((1-D)^2*Ro*C) 
%denominator = (j*omega) ^2 + 1/(Ro*C)*j*omega + ((1-D)^2)/L*C

%Gvd = numerator/denominator
Gvd = tf([0 -Vi/(Ro*C*(1-D)) Vi/(L*C)], [1 1/(Ro*C) (1-D)^2/(L*C)])
bode(Gvd)

grid on

% 360 a  0-
% considerando tudo positivo

%Vo_d = numerator


%G_mod_db = 
% bode = modulo e fase


% frequencia comutacao -> 50 kHz
% analisamos ate a metade da freq de comutacao
% 0 semiplano do tempo ->



% PI só consegue gerar fase de menos 90
% PI tem um polo, entao apos 1e4 nao consigo usar
% 1 POLO GERA FASE DE  -90
% PI = K (s+ w_s) / s