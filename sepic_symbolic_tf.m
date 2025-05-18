%% symbolic test variable


syms Vi Vo D d IL1 IL2 s L1 L2 C1 C2 R;
syms vo_hat d_hat i_L1_hat i_L2_hat v_C1_hat v_C2_hat;

% SEPIC voltage transfer function
Vo_d_relation = Vo == Vi*d/(1-d);
fprintf('SEPIC voltage relationship: ');
disp(Vo_d_relation);

% r output voltage in terms of input voltage and duty cycle
Vo_eq = solve(Vo_d_relation, Vo);
fprintf('Output voltage equation: Vo = ');
disp(Vo_eq);

i_d_hat = (IL1 + IL2)*d_hat + (1-d)*(i_L1_hat + i_L2_hat) ;
v_s_hat = - (Vo)*d_hat + (1-d)*(v_C1_hat + vo_hat) ;

fprintf('Controlled voltage source: v_s_hat = ');
disp(v_s_hat);
fprintf('Controlled current source: i_d_hat = ');
disp(i_d_hat);

% Equations of the SEPIC converter
eq1 = s * L1 * i_L1_hat + (1 - D) * v_C1_hat + (1 - D) * vo_hat - d_hat * Vi - d_hat * Vo == 0;
eq2 = v_C1_hat - s * L2 * i_L2_hat - (1 - D) * v_C1_hat - (1 - D) * vo_hat + d_hat * 2* Vi + d_hat * Vo == 0;
eq3 = s * C1 * v_C1_hat + i_L1_hat - (1 - D) * i_L1_hat - (1 - D) * i_L2_hat + d_hat * (IL1 + IL2) == 0;
eq4 = s * C2 * vo_hat + vo_hat / R - (1 - D) * i_L1_hat - (1 - D) * i_L2_hat + d_hat * (IL1 + IL2) == 0;

% Solve the system
solutions = solve([eq1, eq2, eq3, eq4], [i_L1_hat, i_L2_hat, v_C1_hat, vo_hat]);

% Extract vo_hat / d_hat
if isfield(solutions, 'vo_hat') % Check if there is a solution for vo_hat
    G_symbolic = collect(simplify(solutions.vo_hat / d_hat), s);
end


%%

% Define symbolic variables
syms C1 C2 L1 L2 R IL1 IL2 Vi Vo D s

% Define the numerator and denominator of the transfer function
numerator = [- C1*IL1*L1*L2*R - C1*IL2*L1*L2*R, C1*L1*R*Vi + C1*L2*R*Vi + C1*L1*R*Vo + C1*L2*R*Vo, - C1*L1*R*Vi*D - C1*L2*R*Vi*D - C1*L1*R*Vo*D - C1*L2*R*Vo*D, IL1*L2*R + IL2*L2*R - IL1*L2*R*D - IL2*L2*R*D, - R*Vi - R*Vo + R*Vi*D + R*Vo*D];

denominator = [C1*C2*L1*L2*R, C1*L1*L2, ...
              C1*L1*R + C1*L2*R + C1*L1*R*D^2 + C1*L2*R*D^2, ...
               C2*L1*R*D^2 + C2*L2*R*D^2 - 2*C1*L1*R*D, ...
               -2*C1*L2*R*D - C2*L1*R*D - C2*L2*R*D, ...
              (L1*D^2 - L2*D - L1*D + L2*D^2)*s - R + 2*R*D - R*D^2];

% Substitute numerical values for the parameters
L1_val = 100e-6; % Example: 100 µH
L2_val = 100e-6; % Example: 100 µH
R_val = 10;      % Example: 10 Ω
IL1_val = 1;     % Example: 1 A
IL2_val = 1;     % Example: 1 A
Vi_val = 50  ; % Input voltage [V]
Vo_val = 100 ; % Output voltage [Vo]
D_val = Vo_val / (Vi_val + Vo_val); % SEPIC duty cycle formula
fs = 50 ; % 1/Ts [kHz] 
Ts = (1/fs) * 10-6 ; % secs
Vs  = Vi_val*D_val/(1-D_val);
Po = 250 ; % W
Vrms = sqrt(Vi);
Vpeak = Vrms * sqrt(2);
Amplitude = Vpeak;
Iin = Po/Vi_val ;
IL1_val = Iin;
Io = Po/Vo ;
IL2_val = Io;
Rload = Vo/Io;       % load resistance (Ohm)
R = Rload ;
delta_I_L_p = 0.2 ; % 20 % percentage
delta_I_l1 = I_l1 * delta_I_L_p ;
delta_I_l2 = I_l2 * delta_I_L_p ;
L2 = Vi*D/delta_I_l2*fs ;
L1 = Vi*D/delta_I_l1*fs ;
delta_V_C1_p = .20 ; % 20 % percentage
delta_V_C1 = Vi * delta_V_C1_p; 
delta_Vo = 0.5 ; % 5% 
delta_V_C2 = Vo * delta_Vo; 
C1_val = Io * D / delta_V_C1 * fs;       % Capacitor C1 (F)
C2_val = Io * D / delta_V_C2 * fs;       % Capacitor C1 (F)


% Substitute into expressions
numerator_num = subs(numerator, {C1, C2, L1, L2, R, IL1, IL2, Vi, Vo, D}, {C1_val, C2_val, L1_val, L2_val, R_val, IL1_val, IL2_val, Vi_val, Vo_val, D_val});
denominator_num = subs(denominator, {C1, C2, L1, L2, R, IL1, IL2, Vi, Vo, D}, {C1_val, C2_val, L1_val, L2_val, R_val, IL1_val, IL2_val, Vi_val, Vo_val, D_val});

% Convert to polynomial coefficients
num_coeffs = sym2poly(numerator_num);
den_coeffs = sym2poly(denominator_num);

% Define the transfer function
G = tf(num_coeffs, den_coeffs);

% Display the transfer function
disp(G);

%%

% Kirchhoff's voltage and current law equations in s-domain
eq1 = v_i_hat == s*L1*i_L1_hat + v_C1_hat + (1-d)*v_s_hat ;
eq2 = v_C1_hat == s*L2*i_L2_hat - v_s_hat - s*L1*i_L1_hat ;
eq3 = i_L1_hat == 1/s*L1 * v_s_hat;
eq4 = i_d_hat == vo_hat/R + s*C2*vo_hat;


%iL2_hat = (s*L1*IL1 + Vs_hat + VC1_hat)/s*L2 ;
%iL1_hat = s*C1*VC1_hat + is_hat ;

% Define impedance equivalent (optional approach)
Z_eq = R/(R*s*C2 + 1);
eq4_Zeq = i_d_hat == vo_hat/Z_eq;

fprintf('Using equivalent impedance Z_eqs = %s\n', char(Z_eq));

% System of equations
equations = [eq1, eq2, eq3, eq4_Zeq];
variables = [i_L1_hat, i_L2_hat, v_C1_hat, v_i_hat];

fprintf('Solving system of equations to find vo_hat/d_hat...\n');

% Isolate to find transfer function vo_hat/d_hat (control-to-output)
try
    % Set v_i_hat = 0 to find control-to-output transfer function
    eqs_control = subs(equations, v_i_hat, 0);
    vars_control = [i_L1_hat, i_L2_hat, v_C1_hat];
    sol_control = solve(eqs_control, vars_control);
    
    % Substitute solutions to get vo_hat
    v_o_d_transfer = simplify(vo_hat/d_hat);

    fprintf('Control-to-output transfer function vo_hat/d_hat:\n');
    pretty(v_o_d_transfer);
catch ME
    fprintf('Analytical solution too complex. Using simplified approach.\n');
    
    % Simplified model for SEPIC control-to-output transfer function
    % This is a typical approximation for the SEPIC converter
    simplified_transfer = Vi*(1+s*C1*R)/(1-d)^2;
    fprintf('Simplified control-to-output transfer function (approximate):\n');
    pretty(simplified_transfer);
end

%% Numerical example
fprintf('\nNumerical Example:\n');
fprintf('For Vi = 12V and desired Vo = 24V:\n');


BB = simplify(BB1*d*BB2*(1-d));

% DC MODEL - operation point behaviour
X = - simplify(inv(AA)*BB)*U ;

% Transfer function

E2 = simplify((BB1-BB2)*U)

E = simplify(E1+E2)

C1 = [0 1];

%colect isola variavel, ex em funcao de s
Gvo = collect(simplify(C1*inv(s*eye(2)-AA)*E), s)



