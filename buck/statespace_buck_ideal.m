%%Buck - Modelo CA sem perdas%%

clear all
close all
clc
format long

%% Constantes (componentes)


L = 33e-6;
Cap = 195e-6;
R = 10;
Vin = 36;
D = 0.53572;
fs = 100e3;

%% Compensador
teta=60;
fc = fs/10;
fz = fc*sqrt((1-sind(teta))/(1+sind(teta)));
fp = fc*sqrt((1+sind(teta))/(1-sind(teta)));
fl=fc/10;
%% Variaveis de estado

K = [L 0
    0 Cap];

% Etapa 1
A1 = [0 -1
      1 -1/R];
B1 = [1
      0];
C1 = [1 0 ; 0 1];
E1 = [ 0 ; 0]; 

% Etapa 2
A2 = [0 -1
      1 -1/R];
B2 = [0
      0];
C2 = [1 0 ; 0 1];
E2 = [ 0 ; 0];

%% Definicao das Matrizes Medias

A = D*A1+(1-D)*A2;
B = D*B1+(1-D)*B2;
C = D*C1+(1-D)*C2;
E = D*E1+(1-D)*E2;


% Calculo dos Valores Medios - Considerando o Sistema em Equilibrio

U = Vin;
X = -inv(A)*B*U;
Y = (-C*inv(A)*B+E)*U;

% Veerificacao do ganho estatico


%% Forma Padrao de Representacao em Espaco de Estados
% Nesta, o vetor de estados é X=[i1 i2 v1 v2]' e a entrada U=[u d]'
% O vetor de saida e Y=[i1 i2 v1 v2]'

Ap = inv(K)*A;
Bp = [inv(K)*B inv(K)*((A1-A2)*X+(B1-B2)*U)];
Cp = C;
Ep = [E ((C1-C2)*X+(E1-E2)*U)];

sys=ss(Ap,Bp,Cp,Ep);

%% Conversao para funcao transferencia

[num,den] = ss2tf(Ap,Bp,Cp,Ep,1);
il_vin = tf(num(1,:),den) %il em função de vin
vo_vin = tf(num(2,:),den) %vo em função de vin

[num,den] = ss2tf(Ap,Bp,Cp,Ep,2);
il_d = tf(num(1,:),den) %il em função do duty cycle
vo_d = tf(num(2,:),den) %vo em função do duty cycle - usar este para controle de tensão

fprintf('Vo = %f\n',X(2));
fprintf('Il = %f\n',X(1));
fprintf('fz = %.3fkHz\n',fz*1e-3);
fprintf('fp = %.3fkHz\n',fp*1e-3);
fprintf('fl = %.3fHz\n',fl);
