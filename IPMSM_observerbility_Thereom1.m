% Derivation of the observability matrix for IPMSM
% Reference frame: alpha-beta axes

clear all
close all

syms x1 x2 ua ub x4 x3 Ke R Ld Lq

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
% x3 = rotor angle, x4 = rotor speed
f = A*u - R*A*x - x4*B*x - C*[0; 1]*Ke*x4;
f(3) = x4;
f(4) = 0;

% Output equation h = [ia; ib]
h(1) = x1;
h(2) = x2;

% Initialize observability matrix with 0th Lie derivative
Lf_h = h;
O = jacobian(Lf_h, [x1, x2, x3, x4]);

% Compute higher-order Lie derivatives and extend the observability matrix
for k = 1:3
    % Compute the kth Lie derivative
    Lf_h = jacobian(Lf_h, [x1, x2, x3, x4]) * f;
    % Append its Jacobian to the observability matrix
    O = [O; jacobian(Lf_h, [x1, x2, x3, x4])];
end

syms ud uq id iq 

% Define transformation expressions between alpha-beta and dq frames
ud_expr = ua*cos(x3) + ub*sin(x3);    % ud expression
uq_expr = -ua*sin(x3) + ub*cos(x3);   % uq expression
id_expr = x1*cos(x3) + x2*sin(x3);    % id expression
iq_expr = -x1*sin(x3) + x2*cos(x3);   % iq expression

% Solve substitution rules: express ua, ub, x1, x2 in terms of dq variables
sol = solve([ud == ud_expr, uq == uq_expr, id == id_expr, iq == iq_expr], [ua, ub, x1, x2]);

% Extract substitution expressions
ua_sub = sol.ua;
ub_sub = sol.ub;
x1_sub = sol.x1;
x2_sub = sol.x2;

% Substitute alpha-beta variables with dq variables
O = subs(O, [ua, ub, x1, x2], [ua_sub, ub_sub, x1_sub, x2_sub]);

% Substitute ud, uq with their differential equations
syms did diq
ud_eq = Ld * did + R * id - x4 * Lq * iq;            % ud equation
uq_eq = Lq * diq + R * iq + x4 * Ld * id + Ke * x4;  % uq equation
O = subs(O, [ud, uq], [ud_eq, uq_eq]);

% Simplify using trigonometric identity cos^2 + sin^2 = 1
assume(cos(x3)^2 + sin(x3)^2 == 1);

% Take reduced observability matrix (remove first two rows and columns)
O = O(3:end, 3:end);
O = simplify(O);

% Construct O^T * O
OTO = O.' * O;
OTO = simplify(OTO);

% Determinant of O^T * O
det_OTO = det(OTO);
det_OTO = simplify(det_OTO);

% Polynomial Ds1 (quartic, non-homogeneous), corresponds to D_{s1σ1} in Ref. [1]
Ds1 = subs(det_OTO, [x4], [0]);  

% Verification step: substitute did = diq = 0
% This proves Theorem 1 in Ref. [1]
subs(Ds1, [did, diq], [0, 0])
