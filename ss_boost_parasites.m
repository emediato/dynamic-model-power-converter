clc
clear all
%% eq 1a e 2a etapa
syms C L Ro Vo il vc dvc dil rc rL rs rd vs Vf Vi

% eqn 1a etapa
eqn1 = L*dil == Vi-il*(rc+rs) ;
eqn2 = C*dvc == -Vo/Ro ;

[s1_dil, s1_dvc] = solve(eqn1, eqn2, dil, dvc) ;

% eqn 2a etapa
eqn3 = L*dil == Vi-Vf - il * (rL+rd) - Vo ; 
eqn4 = C*dvc == il - Vo/Ro ; 

[s2_dil, s2_dvc] = solve(eqn3, eqn4, dil, dvc) ; 

s1_dl_isolated = collect(simplify(s1_dil), il);
s1_dvc_isolated = collect(simplify(s1_dvc), vc); % // 0 

s2_dl_isolated = collect(simplify(s2_dil), [vc Vi Vf]);
s2_dvc_isolated = collect(simplify(s2_dvc), vc); % // 0 


%%
% jacobiana das funcoes para pegar o resultado
syms d s 
% dx/dt = Ax + Bu
% [dil/dt ; dvo/dt] = A[il; Vo] + B*U
A1 = [ (-(rc + rs)/L) 0 ; 0 -Ro/(C*Ro) ] ;
A2 = [ Ro*il/(C*Ro) 0 ; 0 -(Vo - Ro*il)/(C*Ro)]
%A1 = [-(rc+rs)/L 0 ; 0 -1/(C*(Ro+rc)) ] ;
%A2 = [-(Ro*rL+Ro*rc)/L 0 ; 0 -1/(C*(Ro+rc)) ] ;

AA = simplify(A1*d+A2*(1-d));
A = AA;

%% matriz B
U = [Vi; Vf] ;

B1 = [1/L 0; 0 0];
B2 = [-(-Ro-rc)/(L*(Ro+rc)) 0 ; -1/L 0];

B = simplify(B1*d+B2*(1-d));

%% MODELO CC
X=-simplify(inv(A)*B)*U 

%FT
E1 = simplify((A1-A2)*X);
E2 = simplify((B1-B2)*U);
E = simplify(E1+E2);
C1 = [0 1];
%%
% Gvd = collect(simplify(C1*inv(s*eye(2)-A )*E,s))

G = C1 / (s * eye(2) - A) * ((A1 - A2) * X + (B1 - B2) * U);




