Vin = 200  ; % Input voltage [V]
Vopk = 311 ; % Output voltage [Vo]

fs = 50e3 ; % 1/Ts [kHz] 
Ts = (1/fs) ; % secs

Dn = 1/2 * (sqrt(4*Vin^2+Vopk^2) - 2*Vin)/(2*Vopk) ; 

dn = Dn;
Vs  = (2*Dn-1)/(Dn*(1-Dn));

Po = 250 ; % W

Vrms = sqrt(Vin);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;

Iin = Po/Vin ;
Io = Po/Vopk ;


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

Gn = subs(G, vars, vals)
%Gn = subs(Gvod, vars, vals)

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
