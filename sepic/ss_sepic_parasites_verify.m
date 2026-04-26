
clc
clear all
syms vc1 vc2 vco dvc1 dvo s
syms C1 C2 L1 L2 Ro Vo il1 il2 Vc dvc dil1 dil2 rc1 rL rL1 rL2 rs rd vs Vf Vi RC rc d ro

%% 
% 1a etapa
Vo=vc2+rc*C2*dvc;

eqn1 = L1*dil1 == Vi-il1*(rL+rs)+rs*il2 ;
eqn2 = L2*dil2 == vc1-rL*il2+rc*il2 ;
eqn3 = C1*dvc1==-il2 ;
eqn4 = C2*dvo == - Vo/(Ro) ; 
[s1_dil1, s1_dil2, s1_dvc1, s1_dvo] = solve(eqn1, eqn2, eqn3, eqn4, dil1, dil2, dvc1, dvo);

s1_dil1_isolated = collect(simplify(s1_dil1), [il1 il2] )
s1_dil2_isolated = collect(simplify(s1_dil2), il2)
s1_dvc1_isolated = collect(simplify(s1_dvc1), vc1) % // 0 
s1_dvo_isolated = collect(simplify(s1_dvo), vc2) % // 0 
%%
% 2a etapa
%io = vc2+ic2+((vc2+ic2*rc2)/Ro) ; 
ic2 = -Vo/Ro +il1+il2 ;
Vl2=-Vf - vc2 -(rL*il2) -ic2*rc -rd*il2 +rd*il2 ;
vc2= ((il1 - il2 - ic2)*Ro)/(ic2*rc)

eqn5 = L1*dil1==Vi-il1*(rL+rc) - rc*il2 - vc1 -Vl2 ;
eqn6 = L2*dil2 == -Vf -il1*(rd+rc) -il2*(rL+rd+rc) +(-rc*Vo)/Ro;
eqn7 = C1*dvc1 == il1 ;
eqn8 = C2*dvo == il1 +il2 -Vo/Ro ;

[s2_dil, s2_dil2, s2_dvc1, s2_dvo] = solve(eqn5, eqn6, eqn7, eqn8,  dil1, dil2, dvc1, dvo);
s2_dl1_isolated = collect(simplify(s2_dil), [il1 il2 Vi Vo])
s2_dl2_isolated = collect(simplify(s2_dil2), [il1 il2 Vi Vf])
s2_dvc1_isolated = collect(simplify(s2_dvc1), vc1) % // 0 
s2_dvo_isolated = collect(simplify(s2_dvo), [il1 il2 Vo]) % // 0 

%%
AA1=[(-rL-rs)/L1 rs/L1 0 0;
    0 (-rL+rc)/L2 1/L2 0;
    0 -1/C1 0 0;
    0 -rc/C2 0 -1/C2*(Ro+rc)];
%AA1 = [0 0 0 0;
%    0 0 1/L2 0;
%    1 - 1/C1 0 0;
%    0 0 0 -1/(Ro*C2)];

AA2=[-(rL+rc)/L1, -rL/L1, -1/L1, -1/L1;
    %-(rd+2*rc)/L2, -(rL2+rd+2*rc1)/L2, 0, -2/L2-rc/(Ro*L2);
    -(rd+rc)/L2, (rd+rc+rL)/L2, 0, -Ro/(L2*(Ro+rc));
    1/C1, 0, 0, 0;
    1/C2, 1/C2, 0, -1/C2*(Ro*rc)];
%AA2 = [0 0 -1/L1 -1/L1;
%        0 0 0 -1/L2;
%        1/L1 0 0 0;
%        0 0 0 -1/(Ro*C2)];

AA = simplify(AA1*d+AA2*(1-d));
A=AA;

BB1 = [1/L1 0 0 0 ; 0  0 0 0; 0  0 0 0; 0  0 0 0];
BB2 = [1/L1 -1/L2 0 0 ; 0 1/L2 0 0; 0 0 0 0; 0  0 0 0];
BB = simplify((BB1*d)+(BB2*(1-d)));

B=BB;
U = [Vi; 0; Vf; 0] ;
%% MODELO CC
% regime permanente
X=-simplify(inv(A)*B)*U 

% FT
E1 = simplify((AA1-AA2)*X);
E2 = simplify((BB1-BB2)*U);
E = simplify(E1+E2);

%
Ca = [0 1 0 0];
Cb = [0 0 1 0];

G1 = Ca / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);
G2 = Cb / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);

G = (simplify(G1+G2))
Gvod = collect(simplify(G), s)

%Gvc1 = collect(simplify(Ca*inv(s*eye(4)-A )*B,s))
%Gvc2 = collect(simplify(Cb*inv(s*eye(4)-A )*B,s))
%Gvo = (simplify(Gvc1+Gvc2))
%% nominal values 
Vin = 50  ; % Input voltage [V]
Vo = 100 ; % Output voltage [Vo]

fs = 50e3 ; % 1/Ts [kHz] 
Ts = (1/fs) ; % secs

Dn = Vo/(Vin+Vo) ; dn = Dn;
Vs  = Vin*Dn/(1-Dn);

Po = 250 ; % W

Vrms = sqrt(Vi);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

Iin = Po/Vin ;
I_l1 = Iin;

Io = Po/Vo ;
I_l2 = Io;

Rload = Vo/Io;       % load resistance (Ohm)
R = Rload ;

delta_I_L_p = 0.2 ; % 20 % percentage
delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;
delta_V_C1_p = .20 ; % 20 % percentage
delta_V_C1 = Vin * delta_V_C1_p; 
delta_Vo = 5/100 ; % 5% 

delta_V_C2 = Vo * delta_Vo; 

rLn = 0.1 ; %100 %mΩ
rcn= 0.010 ; % 10 mΩ
rsn = .08 ;% 80 mΩ
rdn = .08 ;% 80 mΩ
VF = 1 ; % V
Ts  = 1 / fs;
% X == [IL1avg; IL2avg; VC1avg; VC2avg]
IL2_val = Io; I_l2 = IL2_val;
Rload = Vo/Io;       % load resistance (Ohm)
Rn = Rload ;

Iin = Po/Vin ;
IL1_val = Iin; I_l1 = IL1_val;
Vrms = sqrt(Vin);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;
delta_V_C2 = Vo * delta_Vo; 
delta_V_C1 = Vin * delta_Vo; 

% dimensions 

L1_val = (Vin*(dn))/(delta_I_l1*fs) ;
L2_val = (Vin*Dn)/(delta_I_l2*fs) ;
C1_val = (Io * Dn) / (delta_V_C1 * fs);     % Capacitor C1 (F)
C2_val = (Io * Dn) / (delta_V_C2 * fs);    % Capacitor C2 (F)

delta_Q = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;


%% FT real values
rLn = 0.1 ; % 100 %mΩ
rcn= 0.010 ; % 10 mΩ
rsn = .080 ;% 80 mΩ
rdn = .08 ;% 80 mΩ
VFn = 1 ; % V
vars = [Ro  L1  L2  C1  C2  Vi  d rL rc rs rd VF ];
vals = [Rn L1_val L2_val C1_val C2_val Vin Dn rLn rcn rsn rdn VFn ];

%Gn = subs(G, vars, vals)
Gn = subs(Gvod, vars, vals)
[Gn_num, Gn_den] = numden(Gn);
tf_num = double(flip(coeffs(Gn_num)));
tf_den = double(flip(coeffs(Gn_den)));
Gntf = tf(tf_num, tf_den)

%%
step(Gntf)
