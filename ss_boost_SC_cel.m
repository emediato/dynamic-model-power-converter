clc
clear all
syms vc1 vc2 vc3 s
syms C L Ro Vo il vc dvc dil il dvc1 dvc2 dvc3 dvc4 rc rL rs rd vs Vf Vi RC rc d


Vo = vc2+rc*C*dvc2+vc1+rc*C*dvc1 ;

%% 1a etapa
%VL = Vi
eqn1 = L*dil == Vi;
%Vc3 + rc*il2 = Vc1+rc*ic1
eqn2 = vc3+rc*C*dvc3 == vc1+rc*C*dvc1;
% ic2 = -io
eqn3 = C*dvc2==-Vo/Ro;
% ic2 = ic3 + ic1
eqn4 = C*dvc2==C*vc3+C*dvc1;

[s1_dil, s1_dvc1, s1_dvc2, s1_dvc3] = solve(eqn1, eqn2, eqn3, eqn4, dil, dvc1, dvc2, dvc3);

% 2a etapa
eqn5 = L*dil==Vi-vc1-rc*C*dvc1;
eqn6 = vc3+rc*C*dvc3==vc2+rc*C*dvc2;
eqn7 = Vo/Ro +C*dvc2 + C*dvc3 == 0;
%ic2 + ic3 + il = ic1
eqn8 = C*dvc2 + C*dvc3 +il == C*dvc1 ;

[s2_dil, s2_dvc1, s2_dvc2, s2_dvc3] = solve(eqn5, eqn6, eqn7, eqn8, dil, dvc1, dvc2, dvc3)

%%
U = [Vi;0;0;0]
AA1 = [0 0 0 0;
      0 0 0 0; 
      0 0 0 0;
      0 0 0 0]
AA2 = [(2*Ro*Vi + 3*Vi*rc - 2*Ro*vc1 - rc*vc1 + rc*vc2 + rc*vc3 - il*rc^2 - 2*Ro*il*rc)/(L*(2*Ro + 3*rc)) 0 0 0;
      0 -(2*vc1 + vc2 + vc3 - 2*Ro*il - il*rc)/(C*(2*Ro + 3*rc)) 0 0; 
      0 0 0 0;
      0 0 0 0]

AA = simplify(AA1*d+AA2*(1-d));

A=AA;

BB1 = [-1/L 0 0 0 ; 0  0 0 0; 0  0 0 0; 0  0 0 0];
BB2 = [1/L 0 0 0 ; 0  0 0 0; 0  0 0 0; 0  0 0 0];
BB = simplify((BB1*d)+(BB2*(1-d)));

B=BB;

%% MODELO CC
% regime permanente

X=-simplify(inv(A)*B)*U 


%% 
Ca = [0 1 0 0];
Cb = [0 0 1 0];

Gvc1 = collect(simplify(Ca*inv(s*eye(4)-A )*B,s))
Gvc2 = collect(simplify(Cb*inv(s*eye(4)-A )*B,s))
Gvo = (simplify(Gvc1+Gvc2))




