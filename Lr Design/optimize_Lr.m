function [x_opt, P_opt, results] = optimize_Lr(w_max, l_max, params)
% optimize_Lr  Minimize total losses of the resonant inductor
%
% Optimization variables: x = [b, w, k_area, a]
%   b      : center post width/diameter  [m]
%   w      : winding window width        [m]
%   k_area : ratio Ac_legs / Ac_center   [-]
%   a      : center post depth           [m]  (= b for round, free for rectangle)
%
%   params.shape  = 'round' | 'rectangle'
%   params.N_min, params.N_max  integer turn range
%   params.T_max  [°C]  maximum allowed temperature
%   params.R_plate, params.Area_plate  thermal interface params
%   params.T_water, params.sig_grease, params.d_grease

N_range = params.N_min : params.N_max;

P_best = Inf;
x_best = [];
N_best = NaN;

for N = N_range
    [x_N, P_N, flag] = optimize_for_N(N, w_max, l_max, params);
    if flag > 0 && P_N < P_best
        P_best = P_N;
        x_best = x_N;
        N_best = N;
    end
end

if isinf(P_best)
    warning('optimize_Lr: no feasible design found.');
end

x_opt   = [x_best, N_best];
P_opt   = P_best;
results = calculate_Lr(x_best(1), x_best(2), x_best(3), N_best, x_best(4), params.shape, params);

% Append temperature to results
results.T_Lr = calc_Lr_Temp( ...
    results.P_total, results.Ac_center, results.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);

end

% ═════════════════════════════════════════════════════════════════════════
% LOCAL FUNCTIONS
% ═════════════════════════════════════════════════════════════════════════

function [x_opt, P_opt, exitflag] = optimize_for_N(N, w_max, l_max, params)
% x = [b, w, k_area, a]
lb = [1e-3,  2*params.s_ct,  0.1,  1e-3  ];
ub = [100e-3, 100e-3,         2.0,  100e-3];
x0 = [5e-3,  5e-3,           0.5,  5e-3  ];

if strcmp(params.shape, 'round')
    Aeq = [-1, 0, 0, 1];
    beq = 0;
else
    Aeq = [];
    beq = [];
end

objective = @(x) get_Ptotal(x, N, params);
nonlcon   = @(x) constraints(x, N, w_max, l_max, params);

opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

[x_opt, P_opt, exitflag] = fmincon(objective, x0, [], [], Aeq, beq, lb, ub, nonlcon, opts);
end

function P = get_Ptotal(x, N, params)
    r = calculate_Lr(x(1), x(2), x(3), N, x(4), params.shape, params);
    if r.is_valid
        P = r.P_total;
    else
        P = 1e10;
    end
end

function [c, ceq] = constraints(x, N, w_max, l_max, params)
    s_ct = params.s_ct;
    r    = calculate_Lr(x(1), x(2), x(3), N, x(4), params.shape, params);

    % Temperature at this design point
    T_Lr = calc_Lr_Temp( ...
        r.P_total, r.Ac_center, r.h_w, ...
        params.R_plate, params.Area_plate, params.T_water, ...
        params.sig_grease, params.d_grease);

    c(1) = r.Bmax_center - params.Bsat;   % center post flux density
    c(2) = r.Bmax_legs   - params.Bsat;   % outer leg flux density
    c(3) = 2*s_ct - x(2);                 % positive winding width
    c(4) = -r.Vc;                         % Vc > 0
    c(5) = -r.Ac_center;                  % Ac_center > 0
    c(6) = -r.Ac_legs;                    % Ac_legs > 0
    c(7) = r.w_tot - w_max;               % total width limit
    c(8) = r.l     - l_max;               % total length limit
    c(9) = T_Lr    - params.T_max;        % temperature limit

    c(~isfinite(c)) = 1e10;
    ceq = [];
end