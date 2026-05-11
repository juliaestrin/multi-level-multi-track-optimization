function [x_opt, P_opt, results] = optimize_VIRT_2LMT(w_max, l_max, params)
% Optimization variables: x = [b, w, k_area]
%   b      : center post diameter  [m]
%   w      : winding window width  [m]
%   k_area : ratio Ac_legs/Ac_center  [-]

lb = [1e-3,        2*params.s_ct,  0.5];   % lower bounds
ub = [2000e-3,     2000e-3,        0.5];   % upper bounds
x0 = [5e-3,        5e-3,           0.5];   % initial guess (0.5 = old behaviour)

objective = @(x) get_Ptotal(x, params);
nonlcon   = @(x) constraints(x, w_max, l_max, params);

opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
[x_opt, P_opt] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, opts);

% Unpack optimum
b_opt      = x_opt(1);
w_opt      = x_opt(2);
k_area_opt = x_opt(3);

% Full results at optimum
if strcmp(params.centerpost_shape, 'round')
    results = calculateVIRT_round_2LMT(b_opt, w_opt, k_area_opt, params);
else
    % square variant would need similar k_area update — placeholder:
    results = calculateVIRT_square(x_opt(1), x_opt(2), x_opt(3), params);
end

results.topology = params.topology;

% Transformer temperature at optimum
results.T_tx = calculate_transformer_temp( ...
    results.P_total, results.Ac_center, results.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);
end

function P = get_Ptotal(x, params)
b      = x(1);
w      = x(2);
k_area = x(3);

if strcmp(params.centerpost_shape, 'round')
    r = calculateVIRT_round_2LMT(b, w, k_area, params);
else
    r = calculateVIRT_square(b, w, k_area, params);   % update square too when ready
end

if r.is_valid
    P = r.P_total;
else
    P = 1e10;
end
end

function [c, ceq] = constraints(x, w_max, l_max, params)
b      = x(1);
w      = x(2);
k_area = x(3);
s_ct   = params.s_ct;

if strcmp(params.centerpost_shape, 'round')
    r = calculateVIRT_round_2LMT(b, w, k_area, params);
else
    r = calculateVIRT_square(b, w, k_area, params);
end

% Derived widths / lengths
w_winding = w - 2*s_ct;
l         = b + 2*w_winding + 2*s_ct;   % a = b for round

% Transformer temperature
T_tx = calculate_transformer_temp( ...
    r.P_total, r.Ac_center, r.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);

% Inequality constraints: c(i) <= 0
c(1) = r.Bmax_center - params.Bsat;   % center post flux density
c(2) = r.Bmax_legs   - params.Bsat;   % outer leg flux density
c(3) = 2*s_ct - w;                    % positive winding width
c(4) = -r.Vc;                         % Vc > 0
c(5) = -r.Ac_center;                  % Ac_center > 0
c(6) = -r.Ac_legs;                    % Ac_legs > 0
c(7) = r.w_tot - w_max;               % total width constraint
c(8) = l - l_max;                     % length constraint
c(9) = T_tx - params.T_max;           % temperature constraint

ceq = [];
end