
clc
clear all
close all
clc
%%
syms vca vcb vco dvca dvcb dvo dila dilb s d D dvc
syms Ca Cb La Lb Ro Vo ila ilb Vi
% ************************************************************************
% 1st DTs
eqn1 = La*dila == Vi ;
eqn2 = Lb*dilb == Vi - vcb ;
eqn3 = Ca*dvca == -Vo/Ro  ; 
eqn4 = Cb*dvcb == (Vo)/Ro + ilb;
eqn5 =  Vo == vca - vcb;

[s1_dila, s1_dilb, s1_dvca, s1_dvb, s1_dvo] = solve(eqn1, eqn2, eqn3, eqn4, eqn5, dila, dilb, dvca, dvcb, Vo);

s1_dil1_isolated = collect(simplify(s1_dila), [ila ilb] );
s1_dil2_isolated = collect(simplify(s1_dilb), ila);
s1_dvc1_isolated = collect(simplify(s1_dvca), vca) ; % // 0 
s1_dvo_isolated = collect(simplify(s1_dvb), vcb) ; % // 0 
% ************************************************************************
% 2nd (1-D)Ts
eqn6 = La*dila == Vi - vca ;
eqn7 = Lb*dilb == Vi;
eqn8 = Ca*dvca == ila - Vo/Ro ;
eqn9 = Cb*dvcb == Vo/(Ro) ;
eqn10 = Vo == vca - vcb;

[s2_dila, s2_dilb, s2_dvca, s2_dvb, s2_dvo] = solve(eqn6, eqn7, eqn8, eqn9, eqn10, dila, dilb, dvca, dvcb, Vo);

s2_dl1_isolated = collect(simplify(s2_dila), [ila ilb Vi]);
s2_dl2_isolated = collect(simplify(s2_dilb), [ila ilb Vi]);
s2_dvc_isolated = collect(simplify(s2_dvca), vca); % // 0 
s2_dvb_isolated = collect(simplify(s2_dvb), [ila ilb vcb]); % // 0 

%
% ************************************************************************
x = [ila; ilb; vca; vcb];
u = [Vi] ; U = u;

dx1 = [s1_dila; s1_dilb; s1_dvca; s1_dvb];
dx2 = [s2_dila; s2_dilb; s2_dvca; s2_dvb];

AA1 = jacobian(dx1, x)
BB1 = jacobian(dx1, u)

AA2 = jacobian(dx2, x)
BB2 = jacobian(dx2, u)

AA = simplify(AA1*d + AA2*(1-d));
A=AA;

BB = simplify((BB1*d)+(BB2*(1-d)));
B=BB;

%% DC MODEL
% steady state
% X = -simplify(inv(A)*B)*U ;
X = -A\B*U; 

% FT
CA = [0 0 1 0];
CB = [0 0 0 1];

%G1 = CA / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);
G1 = simplify(CA * inv(s*eye(4)-A)*((AA1-AA2)*X+(BB1-BB2)*U));

%G2 = CB / (s * eye(4) - A) * ((AA1 - AA2) * X + (BB1 - BB2) * U);
G2 = simplify(CB * inv(s*eye(4)-A)*((AA1-AA2)*X+(BB1-BB2)*U));

G = (simplify(G1-G2))

Gvod = collect(simplify(G), s)

%% nominal values 
% ************************************************************************
% specifications
Vin = 200; % Input voltage [V]
Vorms = 220 ; 

Vopk = 311 ; % Output voltage [Vo]
fs = 50e3 ; % 1/Ts [kHz] 
Ts  = 1 / fs; 
Po = 250 ; % W

Virms = Vin * sqrt(2);
Vopkrms = Vopk * sqrt(2);

Vop = Vorms * sqrt(2) ; 
Iop = 2*Po/Vop ; 

Io = Po/Vopk ;
Iop = (Po / Vop) * (2) ;

Rload = (Vorms^2)/Po; Rn = Rload;  % load resistance (Ohm)

% ************************************************************************
% duty cycle
% when theta = 90
n = 90*pi/180 ;

Dn = 1/2 + (sqrt(4*Vin^2 + (Vop*sin(n))^2) - 2*Vin)/(2*Vop*sin(n))
dn = Dn;
% Dn = 1/2 * (sqrt(4*Vin^2+Vopk^2) - 2*Vin)/(2*Vopk) ; dn = Dn;

dDC = 1/2 ;
dAC = (sqrt(4*Vin^2+Vop^2) - 2*Vin)/(2*Vop) ;

% max duty cycle
Dmax = dDC + dAC ;

% ************************************************************************
% inductors
% delta_iL max occurs in 90 degrees
Ilpk = Iop / (1-Dmax) ;
Ila = Io/(1-Dn) ;

e_iL = 25/100 ;

La = Vin *(Dmax) / (fs*Ilpk*e_iL) ; L1_val = La ;  L2_val  = La ; 
Lb = La;

% ************************************************************************
% capacitors
e_VC = 1/100 ;

Vapk = Vin/(1-Dmax) ;
Ca = (Iop*Dmax) / (fs*Vapk*e_VC) ;

C1_val = Ca ; C2_val = C1_val ; Cb = C2_val ; 

% ************************************************************************
Vo_avg  = (2*Dn-1)/(Dn*(1-Dn))/Vin ;

%% FT real values
vars = [Ro La Lb Ca Cb Vi Vo d D];
vals = [Rn L1_val L2_val C1_val C2_val Vin Vopk Dn Dn];

Gn = subs(G, vars, vals)
% Gn = subs(Gvod, vars, vals)

[Gn_num, Gn_den] = numden(-(105967757824964927728555231629219735061127221305409536*s^3 - 7720861047202367512292233219647958879185247808854276702208*s^2 + 10610893122177765069019208208487502455088642307254580951384064*s - 277167904149427990983735043322923553340547265573034687090013503488)/(150247203216521486179182832697557873729606778880*s^4 + 437883296797375696914773586802206551276680233615360*s^3 + 10787320783327900772600382798024309177123343857654693888*s^2 + 15719386075384384808284907076684636608632396833233420943360*s + 120638453862536964755128727741200042858109671565762559140167680));
tf_num = double(flip(coeffs(Gn_num)));
tf_den = double(flip(coeffs(Gn_den)));
Gntf = tf(tf_num, tf_den)


step(Gntf)
