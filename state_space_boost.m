syms d Vi C L Ro s


% Variáveis de estado
U = [Vi;0];

AA1 = [0 0; 0 -1/(Ro*C)];

AA2 = [0 -1/L; 1/C -1/(Ro*C)];

AA = simplify(AA1*d+AA2*(1-d));

BB1 = [1/L 0; 0 0];

BB2 = [1/L 0; 0 0];

BB = simplify(BB1*d+BB2*(1-d));


%% CC MODEL - operation point behaviour
X = - simplify(inv(AA)*BB)*U ;

%% Transfer function
E1 = ( AA1 - AA2)*X

E2 = simplify((BB1-BB2)*U)

E = simplify(E1+E2)

C1 = [0 1];

%colect isola variavel, ex em funcao de s
Gvo = collect(simplify(C1*inv(s*eye(2)-AA)*E), s)

