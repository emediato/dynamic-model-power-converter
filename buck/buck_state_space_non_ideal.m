%%Buck - Modelo CA com perdas%%

%clear all
%close all
clc
format long

%% Constantes (componentes)

L = 500e-6;
Cap = 10e-6;
R = 35; %
Vin = 36;
D = 0.5;

%Constantes das perdas
Rl=0.135;
Vd=0.6;
Rs=0.07;

%% Variaveis de estado

K = [L 0
    0 Cap];

% Etapa 1
A1 = [(Rs-Rl) -1
      1 -1/R];
B1 = [1 0
     0 0];
C1 = [1 0 ; 0 1];
E1 = [0 0 ; 0 0]; 

% Etapa 2
A2 = [-Rl -1
      1 -1/R];
B2 = [0 -1
      0 0];
C2 = [1 0 ; 0 1];
E2 = [0 0 ; 0 0];

%% Definicao das Matrizes Medias

A = D*A1+(1-D)*A2;
B = D*B1+(1-D)*B2;
C = D*C1+(1-D)*C2;
E = D*E1+(1-D)*E2;


% Calculo dos Valores Medios - Considerando o Sistema em Equilibrio

U = [Vin;Vd];
X = -inv(A)*B*U;
Y = (-C*inv(A)*B+E)*U;



%% Forma Padrao de Representacao em Espaco de Estados


Ap = inv(K)*A;
Bp = [inv(K)*B inv(K)*((A1-A2)*X+(B1-B2)*U)];
Cp = C;
Ep = [E ((C1-C2)*X+(E1-E2)*U)];

sys=ss(Ap,Bp,Cp,Ep);

%% Conversao para funcao transferencia

[num,den] = ss2tf(Ap,Bp,Cp,Ep,1);
il_vin = tf(num(1,:),den)
vo_vin = tf(num(2,:),den)

[num,den] = ss2tf(Ap,Bp,Cp,Ep,2);
il_d = tf(num(1,:),den)
vo_d = tf(num(2,:),den)

fprintf('Vo = %f\n',X(2));
fprintf('Il = %f\n',X(1));
