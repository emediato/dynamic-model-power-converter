
% Definição de variáveis simbólicas
syms dil1 dil2 Dd dvc1 dvc2 s D il2 il1 Vout Vin Vc1 Vc2 C2sym C1sym L2sym L1sym Rsym 

%% Etapa 1 (chave ligada)
% Equações do inversor Boost na primeira etapa
eq1 =  L1sym * dil1 == Vin;
eq2 =  L2sym * dil2 == Vin - Vc2;
eq3 =  C1sym * dvc1 == -Vout/Rsym; 
eq4 =  C2sym * dvc2 == Vout/Rsym + il2;
eq5 =  Vout == Vc1 - Vc2;

% Resolver o sistema
[s1_dil1, s2_dil2, s3_dvc1, s4_dvc2, s5_Vout] = solve(eq1, eq2, eq3, eq4, eq5, dil1, dil2, dvc1, dvc2, Vout);

%% Etapa 2 (chave desligada)
% Equações do inversor Boost na segunda etapa
eq6 =  L1sym * dil1 == Vin - Vc1;
eq7 =  L2sym * dil2 == Vin;
eq8 =  C1sym * dvc1 == -Vout/Rsym + il1;
eq9 =  C2sym * dvc2 == Vout/Rsym;
eq10 =  Vout == Vc1 - Vc2;

% Resolver o sistema
[s6_dil1, s7_dil2, s8_dvc1, s9_dvc2, s10_Vout] = solve(eq6, eq7, eq8, eq9, eq10, dil1, dil2, dvc1, dvc2, Vout);

% Vetor de entradas
U = Vin;

% Vetor de funções da etapa 1
f1 = [s1_dil1; s2_dil2; s3_dvc1; s4_dvc2];

% Vetor de funções da etapa 2
f2 = [s6_dil1; s7_dil2; s8_dvc1; s9_dvc2];

% Vetor de estados
x = [il1; il2; Vc1; Vc2];

% Jacobiano da etapa 1
A1 = jacobian(f1, x);

% Jacobiano da etapa 2
A2 = jacobian(f2, x);

% Matriz B da etapa 1 (derivadas das equações em relação às entradas)
B1 = jacobian(f1, U);

% Matriz B da etapa 2
B2 = jacobian(f2, U);

% Matrizes médias ponderadas
A = A1*D + A2*(1-D);
B = B1*D + B2*(1-D);

%% MODELO CC

X = -A\B*U;  % Troquei inv(A) por A\B

%% MATRIZ C

C1mat = [0 0 1 0];
C2mat = [0 0 0 1];

%% Funções de Transferência

GC1 = simplify(C1mat*inv(s*eye(4)-A)*((A1-A2)*X+(B1-B2)*U));
GC2 = simplify(C2mat*inv(s*eye(4)-A)*((A1-A2)*X+(B1-B2)*U));

Gvout = simplify(GC1 - GC2);

disp('A1 = '); disp(A1);
disp('A2 = '); disp(A2);
disp('B1 = '); disp(B1);
disp('B2 = '); disp(B2);

%% Substituição de valores

% Definir valores
C1num = 22e-6;
C2num = 22e-6;
L1num = 2.2e-3;
L2num = 2.2e-3;
Vinnum = 200;
Dnum = 0.672;
Rnum = 193.4; % Valor do resistor de carga, defina conforme seu caso!

% Substituir nas equações
Gvout_num = subs(Gvout, [C1sym, C2sym, L1sym, L2sym, Vin, D, Rsym], [C1num, C2num, L1num, L2num, Vinnum, Dnum, Rnum]);

% Separar numerador e denominador
[num, den] = numden(Gvout_num);



% Exibir coeficientes em notação de engenharia
disp('Numerador:');
disp(vpa(coeffs(Gvonum, s, 'All'), 4));

disp('Denominador:');
disp(vpa(coeffs(Gvoden, s, 'All'), 4));
