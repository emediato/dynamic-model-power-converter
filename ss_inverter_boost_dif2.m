
clc
clear all
%%
syms vca vcb vco dvca dvcb dvo dila dilb s d D dvc
syms Ca Cb La Lb Ro Vo ila ilb Vi

%
% 1a etapa

eqn1 = La*dila == Vi  ;
eqn2 = Lb*dilb == Vi - vcb ;

eqn3 = Cb*dvcb == -(vca-vcb)/Ro ;
eqn4 = Ca*dvca == (vca-vcb)/Ro + ilb %-ilb  ; 


[s1_dila, s1_dilb, s1_dvca, s1_dvb] = solve(eqn1, eqn2, eqn3, eqn4, dila, dilb, dvca, dvcb);

s1_dil1_isolated = collect(simplify(s1_dila), [ila ilb] )
s1_dil2_isolated = collect(simplify(s1_dilb), ila)
s1_dvc1_isolated = collect(simplify(s1_dvca), vca) % // 0 
s1_dvo_isolated = collect(simplify(s1_dvb), vcb) % // 0 


%
% 2a etapa
eqn5 = La*dila == Vi - vca ;
eqn6 = Lb*dilb == Vi;

eqn7 = Ca*dvca == ila - (vcb)/(Ro) ;
eqn8 = Cb*dvcb == (vcb)/(Ro) ;


[s2_dila, s2_dilb, s2_dvca, s2_dvb] = solve(eqn5, eqn6, eqn7, eqn8, dila, dilb, dvca, dvcb);

s2_dl1_isolated = collect(simplify(s2_dila), [ila ilb Vi])
s2_dl2_isolated = collect(simplify(s2_dilb), [ila ilb Vi])
s2_dvc_isolated = collect(simplify(s2_dvca), vca) % // 0 
s2_dvb_isolated = collect(simplify(s2_dvb), [ila ilb vcb]) % // 0 

%%

x = [ila; ilb; vca; vcb];
u = [Vi] ; U = u;

dx1 = [s1_dila; s1_dilb; s1_dvca; s1_dvb];
dx2 = [s2_dila; s2_dilb; s2_dvca; s2_dvb];

AA1 = jacobian(dx1, x)
BB1 = jacobian(dx1, u)

AA2 = jacobian(dx2, x)
BB2 = jacobian(dx2, u)

AA = simplify(AA1*d+AA2*(1-d));
A=AA;

BB = simplify((BB1*d)+(BB2*(1-d)));

B=BB;

%% MODELO CC
% regime permanente
X=-simplify(inv(A)*B)*U 

% FT
E1 = simplify((AA1-AA2)*X);
E2 = simplify((BB1-BB2)*U);
E = simplify(E1+E2);

%
CA = [0 0 0 1];
%Cb = [0 0 1 0];

G1 = CA / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);
%G = (simplify(G1+G2))

G = (simplify(G1))

Gvod = collect(simplify(G), s)

%Gvc1 = collect(simplify(Ca*inv(s*eye(4)-A )*B,s))
%Gvc2 = collect(simplify(Cb*inv(s*eye(4)-A )*B,s))
%Gvo = (simplify(Gvc1+Gvc2))
%% nominal values 
Vin = 200  ; % Input voltage [V]
Von = 311 ; % Output voltage [Vo]

fs = 50e3 ; % 1/Ts [kHz] 
Ts = (1/fs) ; % secs

Dn = 1/2 * (sqrt(4*Vin^2+Von^2) - 2*Vin)/(2*Von) ; 

dn = Dn;
Vs  = (2*Dn-1)/(Dn*(1-Dn));

Po = 250 ; % W

Vrms = sqrt(Vin);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

Iin = Po/Vin ;
I_l1 = Iin;

Io = Po/Von ;
I_l2 = Io;

Rload = Von/Io;       % load resistance (Ohm)
R = Rload ;

delta_I_L_p = 0.2 ; % 20 0% percentage
delta_V_C1_p = .20 ; % 20 % percentage
delta_Vo = 5/100 ; % 5% 
Ts  = 1 / fs;

% X == [IL1avg; IL2avg; VC1avg; VC2avg]

IL2_val = Io; I_l2 = IL2_val;
Rload = Von/Io;       % load resistance (Ohm)
Rn = Rload ;

Iin = Po/Vin ;
IL1_val = Iin; I_l1 = IL1_val;
Vrms = sqrt(Vin);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;
delta_V_C2 = Von * delta_Vo; 
delta_V_C1 = Vin * delta_Vo; 

% dimensions 
Ila_val = Io/(1-D);
L1_val = (Vin*(dn))/(delta_I_l1*fs) ;
L2_val = (Vin*Dn)/(delta_I_l2*fs) ;
C1_val = (Io * Dn) / (delta_V_C1 * fs);     % Capacitor C1 (F)
C2_val = (Io * Dn) / (delta_V_C2 * fs);    % Capacitor C2 (F)

delta_Q = (Ts/2 * delta_I_L_p/2)/2;
C = delta_Vo/delta_Q ;


%% FT real values
vars = [Ro La Lb Ca Cb Vi Vo d];
vals = [Rn L1_val L2_val C1_val C2_val Vin Von Dn];

Gn = subs(Gvod, vars, vals)
%Gn = subs(G, vars, vals)

[Gn_num, Gn_den] = numden(Gn);
tf_num = double(flip(coeffs(Gn_num)));
tf_den = double(flip(coeffs(Gn_den)));
Gntf = tf(tf_num, tf_den)


%%
[Gvonum,Gvoden] = numden(Gn);
Gvonumcoefs = flip(double(coeffs(Gvonum)));
Gvodencoefs = flip(double(coeffs(Gvoden)));

rscl = Gvodencoefs(end);
Gvonumcoefs = Gvonumcoefs / rscl;
Gvodencoefs = Gvodencoefs / rscl;

gain = Gvonumcoefs(end);
Gvonumcoefs = Gvonumcoefs / gain;

format shortG;

fprintf('gain: %.3g\nnum: ', gain);
fprintf(' %.3g', [zeros(1, length(Gvodencoefs)-length(Gvonumcoefs)) Gvonumcoefs]);
fprintf('\nden: ');
fprintf(' %.3g', Gvodencoefs);
fprintf('\n');


%%
step(Gntf)