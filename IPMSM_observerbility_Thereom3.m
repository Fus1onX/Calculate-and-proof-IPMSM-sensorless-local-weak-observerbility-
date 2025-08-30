% Derivation of the observability matrix for the nonlinear IPMSM model
% Reference frame: alpha-beta axes

clear all
close all

syms x1 x2 ua ub x3 we Ke R Ld Lq

% Matrix A
A = (Ld - Lq)/(2*Ld*Lq)*[-cos(2*x3), -sin(2*x3); 
                          -sin(2*x3), cos(2*x3)] + ...
    (Ld + Lq)/(2*Ld*Lq)*eye(2);

% Matrix B
B = (Ld^2 - Lq^2)/(2*Ld*Lq)*[-sin(2*x3), cos(2*x3); 
                               cos(2*x3), sin(2*x3)] + ...
    ((Ld - Lq)^2)/(2*Ld*Lq)*[0, -1; 1, 0];

% Matrix C
C = [cos(x3), -sin(x3); 
     sin(x3),  cos(x3)] * ...
    [1/Ld, 0; 0, 1/Lq];

% State vector x = [ia; ib]
x = [x1; x2];

% Input vector u = [u_alpha; u_beta]
u = [ua; ub];

% State-space model: dx/dt = f(x,u)
% x3 = rotor angle, we = electrical speed
% Note: substituting we = 0 does not change det(OTO)
f = A*u - R*A*x - we*B*x - [-sin(x3)/Lq ; cos(x3)/Lq]*Ke*we;
f(3) = we;

% Output equation h = [ia; ib]
h(1) = x1;
h(2) = x2;

% Initialize observability matrix with 0th Lie derivative
Lf_h = h;
O = jacobian(Lf_h, [x1, x2, x3]);

% Compute Lie derivatives up to order 2 and expand the observability matrix
for k = 1:2
    % kth Lie derivative
    Lf_h = jacobian(Lf_h, [x1, x2, x3]) * f;
    % Append its Jacobian to O
    O = [O; jacobian(Lf_h, [x1, x2, x3])];
end

syms ud uq id iq 

% Define transformations between alpha-beta and dq frames
ud_expr = ua*cos(x3) + ub*sin(x3);
uq_expr = -ua*sin(x3) + ub*cos(x3);
id_expr = x1*cos(x3) + x2*sin(x3);
iq_expr = -x1*sin(x3) + x2*cos(x3);

% Solve substitution rules: express ua, ub, x1, x2 in terms of dq variables
sol = solve([ud == ud_expr, uq == uq_expr, id == id_expr, iq == iq_expr], [ua, ub, x1, x2]);

% Extract substitution expressions
ua_sub = sol.ua;
ub_sub = sol.ub;
x1_sub = sol.x1;
x2_sub = sol.x2;

% Apply substitution: replace alpha-beta variables with dq variables
O = subs(O, [ua, ub, x1, x2], [ua_sub, ub_sub, x1_sub, x2_sub]);

% Simplify using trigonometric identity cos^2 + sin^2 = 1
assume(cos(x3)^2 + sin(x3)^2 == 1);
O = simplify(O);

% Substitute ud, uq with their dynamic equations
syms did diq
ud_eq = Ld * did + R * id - we * Lq * iq;
uq_eq = Lq * diq + R * iq + we * Ld * id + Ke * we;
O = subs(O, [ud, uq], [ud_eq, uq_eq]);

% Set we = 0 for analysis
O = subs(O, [we], [0]);

% Simplify observability matrix
assume(cos(x3)^2 + sin(x3)^2 == 1);
O = simplify(O);

% Construct O^T * O
OTO = O.' * O;

% Determinant of O^T * O
% This corresponds to D_sigma2 in the paper (Theorem 3)
D_sigma2 = simplify(det(OTO));
