% Derivation of the observability matrix for IPMSM nonlinear equations
% Simplified from 8×4 to 6×2 form to reduce expression complexity
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

% Compute higher-order Lie derivatives and expand observability matrix
for k = 1:3
    % kth Lie derivative
    Lf_h = jacobian(Lf_h, [x1, x2, x3, x4]) * f;
    % Append Jacobian of the Lie derivative
    O = [O; jacobian(Lf_h, [x1, x2, x3, x4])];
end

syms ud uq id iq 

% Transformation between alpha-beta and dq frames
ud_expr = ua*cos(x3) + ub*sin(x3);
uq_expr = -ua*sin(x3) + ub*cos(x3);
id_expr = x1*cos(x3) + x2*sin(x3);
iq_expr = -x1*sin(x3) + x2*cos(x3);

% Solve for substitution rules: express ua, ub, x1, x2 in terms of dq variables
sol = solve([ud == ud_expr, uq == uq_expr, id == id_expr, iq == iq_expr], [ua, ub, x1, x2]);

% Extract substitution expressions
ua_sub = sol.ua;
ub_sub = sol.ub;
x1_sub = sol.x1;
x2_sub = sol.x2;

% Apply substitution: replace alpha-beta with dq variables
O = subs(O, [ua, ub, x1, x2], [ua_sub, ub_sub, x1_sub, x2_sub]);

% Substitute ud, uq with electrical equations
syms did diq
ud_eq = Ld * did + R * id - x4 * Lq * iq;
uq_eq = Lq * diq + R * iq + x4 * Ld * id + Ke * x4;
O = subs(O, [ud, uq], [ud_eq, uq_eq]);

% Simplify with trigonometric identity cos^2 + sin^2 = 1
assume(cos(x3)^2 + sin(x3)^2 == 1);

% Reduced observability matrix: remove first two rows and columns
O = O(3:end, 3:end);
O = simplify(O);

% Construct O^T * O
OTO = O.' * O;
OTO = simplify(OTO);

% Determinant of O^T * O
det_OTO = det(OTO);

% Polynomial Ds1 (quartic, non-homogeneous), corresponds to D_{s1σ1} in Ref. [1]
Ds1 = subs(det_OTO, [x4], [0]);
Ds1 = simplify(Ds1);

% Assume positive motor parameters
assumeAlso([Ke Ld Lq R] > 0)

% Only numerator is considered
Dss1 = Ds1*(Ld^10*Lq^10);

% --- Extract quadratic form coefficients: A, b, c ---
a11 = 0.5*diff(Dss1, id, 2);
a22 = 0.5*diff(Dss1, iq, 2);
a12 = 0.5*diff(diff(Dss1, id), iq);
G   = [a11, a12; a12, a22];

b1  = subs(diff(Dss1, id), [id,iq], [0,0]);
b2  = subs(diff(Dss1, iq), [id,iq], [0,0]);
b   = [b1; b2];

c   = subs(Dss1, [id,iq], [0,0]);

G = simplify(G);

% General expression (using pseudoinverse)
kappa = simplify( c - sym(1)/4 * ( b.' * pinv(G) * b ) );  
% This verifies the expressions for J1, J2, J3 in the referenced paper

% --------------------------------------------------------
% Manual extraction of J2 part for further analysis
% Known: J2 = A*x^2 + B*x*y + C*y^2, where x = did^2, y = diq^2
% --------------------------------------------------------
assumeAlso([Ld Lq R] > 0)
x = did^2;
y = diq^2;

% J2 polynomial (manually extracted)
J2 = (4*Ld^6*Lq^4*x*y + 8*Ld^6*Lq^2*R^2*x*y + 4*Ld^6*R^4*x*y ...
    - 8*Ld^5*Lq^5*x*y + 8*Ld^5*Lq*R^4*x*y + 4*Ld^4*Lq^6*x*y ...
    + Ld^4*Lq^4*R^2*x^2 - 14*Ld^4*Lq^4*R^2*x*y + Ld^4*Lq^4*R^2*y^2 ...
    + 9*Ld^4*Lq^2*R^4*x*y + Ld^4*Lq^2*R^4*y^2 + 13*Ld^4*R^6*x*y + Ld^4*R^6*y^2 ...
    - 40*Ld^3*Lq^3*R^4*x*y + 8*Ld^2*Lq^6*R^2*x*y + Ld^2*Lq^4*R^4*x^2 ...
    + 9*Ld^2*Lq^4*R^4*x*y - 24*Ld^2*Lq^2*R^6*x*y + 9*Ld^2*R^8*x*y ...
    + 8*Ld*Lq^5*R^4*x*y - 18*Ld*Lq*R^8*x*y + 4*Lq^6*R^4*x*y ...
    + Lq^4*R^6*x^2 + 13*Lq^4*R^6*x*y + 9*Lq^2*R^8*x*y);

syms x y real
J2 = simplify( subs(J2, [did^2, diq^2], [x, y]) );
J2 = collect(J2, [x y]);   % Organize as polynomial in x, y

% --- Extract A, B, C (quadratic form coefficients) ---
j1 = simplify( diff(J2, x, 2) / 2 );   % coefficient of x^2
j2 = simplify( diff(J2, y, 2) / 2 );   % coefficient of y^2
j3 = simplify( diff(diff(J2, x), y) ); % coefficient of x*y

% % Optional: check consistency
% J2_check = simplify( J2 - (j1*x^2 + j2*x*y + j3*y^2) );  % should be 0

% --- Key expression: B + 2*sqrt(A*C) ---
expr = simplify( j3 ); % (optionally + 2*sqrt(A*C))

% --------------------------------------------------------
% Numerical verification: ensure j3 > 0 over parameter ranges
% --------------------------------------------------------

% 1) Create function handle
f_expr = matlabFunction(expr, 'Vars', [Ld, Lq, R]);

% 2) Define parameter ranges
Ld_vals = 1e-3:0.5e-3:30e-3;  % Ld [H]
Lq_vals = 1e-3:0.5e-3:30e-3;  % Lq [H]
R_vals  = 0.01:0.1:10;        % R [Ohm]

% 3) Generate 3D mesh grid
[LD, LQ, RR] = ndgrid(Ld_vals, Lq_vals, R_vals);

% 4) Vectorized evaluation
try
    V = f_expr(LD, LQ, RR);
catch
    % If function only supports scalar input, fall back to arrayfun
    V = arrayfun(f_expr, LD, LQ, RR);
end

% 5) Positivity check
isPositive = all(V(:) > 0);
if ~isPositive
    idx = find(V <= 0, 1, 'first');
    [i1, i2, i3] = ind2sub(size(V), idx);
    fprintf('Non-positive point found: Ld=%g, Lq=%g, R=%g, val=%g\n', ...
        Ld_vals(i1), Lq_vals(i2), R_vals(i3), V(i1,i2,i3));
else
    fprintf('All sampled points are positive.\n');
end

if isPositive
    disp('✅ Within sampled ranges, expr > 0 always holds.');
else
    disp('❌ Found point(s) where expr <= 0.');
end
