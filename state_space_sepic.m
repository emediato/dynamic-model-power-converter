syms d Vi C1 C2 L1 L2 Ro s


% Variáveis de estado
U = [Vi;0;0;0];

AA1 = [0 0 0 0;
    0 0 1/L2 0;
    1 - 1/C1 0 0;
    0 0 0 -1/(Ro*C2)];

AA2 = [0 0 -1/L1 -1/L1;
        0 0 0 -1/L2;
        1/L1 0 0 0;
        0 0 0 -1/(Ro*C2)];

AA = simplify(AA1*d+AA2*(1-d));
A = AA;

BB1 = [1/L1 0; 0 0];

BB2 = [1/L1 0; 0 0];

BB = simplify(BB1*d*BB2*(1-d));

X = -simplify(BB)