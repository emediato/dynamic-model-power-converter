
clc
clear all
%%
syms VA va VB vb VC Vc dvca dvcb dvco dila dilb s d D dvc
syms Ca Cb La Lb Ro Vo Vi ILB ILA ila ilb io

% Modelo de pequenos sinais: descreve a dinâmica do conversor
eqn1 = ila*(1-D) - ILA*d == va*s*Ca + io  ;
eqn2 = s * La * ila  == -va*(1-D) + va*d ;
eqn3 = ilb * d + ilb *(D) + io == - s * Cb * vb ; 
eqn4 =  s * Lb * ilb  == - VB*d -vb*D ;

[s1_dila, s1_dilb, s1_dvca, s1_dvcb] = solve(eqn1, eqn2, eqn3, eqn4, va, ila, vb, ilb);

s1_dil1_isolated = collect(simplify(s1_dila), [ila] )
s1_dil2_isolated = collect(simplify(s1_dilb), ilb)
s1_dvc1_isolated = collect(simplify(s1_dvca), va) % // 0 
s1_dvo_isolated = collect(simplify(s1_dvcb), vb) % // 0 


%%
x = [ila; ilb; va; vb];
u = [Vi] ; U = u;
