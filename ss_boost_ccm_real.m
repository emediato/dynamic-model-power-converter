clc
clear all
close all

%% Equações 1ª e 2ª etapas
%% 1ª etapa
syms dvc dil vc il Vi L C Ro rl rd rc rs VF 
Vo=vc+rc*C*dvc;
eqn1=L*dil==Vi-il*(rs+rl);
eqn2=C*dvc==-Vo/Ro;

[s1_dil,s1_dvc]=solve(eqn1,eqn2,dil,dvc)

%% 2ª etapa
eqn3=L*dil==Vi-Vo-il*(rd+rl)-VF;
eqn4=C*dvc==il-Vo/Ro;

[s2_dil,s2_dvc]=solve(eqn3,eqn4,dil,dvc)


%% Variáveis de entrada
syms vc il Vi L C Ro rl rd rc rs VF d s

U=[Vi; VF];

%% Matriz A
AA1=[-(rl+rs)/L        0;
      0          -1/(C*(Ro+rc))];
    
AA2=[ -(Ro*rc+Ro*rd+Ro*rl+rc*rd+rc*rl)/(L*(Ro+rc))   -Ro/(L*(Ro+rc));
                     Ro/(C*(Ro+rc))                   -1/(C*(Ro+rc))];
       

AA=simplify(AA1*d + AA2*(1-d));
 
%% Matriz B
BB1=[  1/L  0;
        0    0];
BB2=[  -(-Ro-rc)/(L*(Ro+rc))   -1/L;
         0                       0];

BB=simplify(BB1*d + BB2*(1-d));

%% Modelo CC - Descreve Conversor no Ponto de Operação

A = AA;
X = -simplify(inv(AA)*BB)*U


%% Função de Transferência

E1 = simplify((AA1-AA2)*X);
E2 = simplify((BB1-BB2)*U);

E = simplify(E1+E2);

%vo/d
C1 = [0 1];  %Coloca 'um' na variavel de estado, na qual, deseja-se a função de transferência
Gvo = collect(simplify(C1*inv(s*eye(2)-A)*E), s)  %eye é a matriz identidade


