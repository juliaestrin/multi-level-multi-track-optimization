function [x_opt, P_opt, results] = optimize_Lr_EQ(w_max, l_max, params)
% optimize_Lr_EQ  Minimize total losses of the EQ resonant inductor.
%
% Optimization variables: x = [b, w, k_area]
%   b      : center post diameter      [m]
%   w      : winding window width      [m]
%   k_area : Ac_center / Ac_legs       [-]
%
% Integer sweep: N in [params.N_min, params.N_max]

N_range = params.N_min : params.N_max;
P_best  = Inf;
x_best  = [];
N_best  = NaN;

for N = N_range
    [x_N, P_N, flag] = optimize_for_N(N, w_max, l_max, params);
    fprintf('N=%d  P=%.3f W  flag=%d\n', N, P_N, flag);
    if flag > 0 && P_N < P_best
        P_best = P_N;
        x_best = x_N;
        N_best = N;
    end
end

if isinf(P_best)
    warning('optimize_Lr_EQ: no feasible design found across N = %d..%d.', ...
            params.N_min, params.N_max);
    x_opt   = [NaN, NaN, NaN, NaN];
    P_opt   = Inf;
    results = struct('is_valid', false, 'fail_reason', 'No feasible design found');
    return;
end

x_opt   = [x_best, N_best];
P_opt   = P_best;
results = calculate_Lr_EQ(x_best(1), x_best(2), x_best(3), N_best, params);
results.N    = N_best;
results.T_Lr = calc_Lr_Temp( ...
    results.P_total, results.Ac_center, results.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);
end


% =========================================================================
function [x_opt, P_opt, exitflag] = optimize_for_N(N, w_max, l_max, params)
% b, w, k
lb = [1e-3,            2*params.s_ct + 0.5e-3,  0.1];
ub = [w_max,         w_max,                  2.0];
x0 = [10e-3, 10e-3,         0.5];

opts = optimoptions('fmincon', ...
    'Algorithm',  'sqp', ...
    'Display',    'off');opts = optimoptions('fmincon', ...
    'Algorithm',           'sqp', ...
    'Display',             'off', ...
    'MaxIterations',       1000,  ...       % up from default 400
    'MaxFunctionEvaluations', 10000, ...    % up from default 3000
    'OptimalityTolerance', 1e-6,  ...
    'StepTolerance',       1e-8);

[x_opt, P_opt, exitflag] = fmincon( ...
    @(x) objective(x, N, params), ...
    x0, [], [], [], [], lb, ub, ...
    @(x) constraints(x, N, w_max, l_max, params), ...
    opts);
end


% =========================================================================
function P = objective(x, N, params)
r = calculate_Lr_EQ(x(1), x(2), x(3), N, params);
if r.is_valid
    P = r.P_total;
else
    P = 1e10;
end
end


% =========================================================================
function [c, ceq] = constraints(x, N, w_max, l_max, params)
r = calculate_Lr_EQ(x(1), x(2), x(3), N, params);

if ~r.is_valid
    c   = 1e10 * ones(9,1);
    ceq = [];
    return;
end

T_Lr = calc_Lr_Temp( ...
    r.P_total, r.Ac_center, r.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);

c(1) = r.Bmax_center - params.Bsat;
c(2) = r.Bmax_legs   - params.Bsat;
c(3) = 2*params.s_ct - x(2);
c(4) = -r.Vc;
c(5) = -r.Ac_center;
c(6) = -r.Ac_legs;
c(7) = r.w_tot - w_max;
c(8) = r.l     - l_max;
c(9) = T_Lr    - params.T_max;

c(~isfinite(c)) = 1e10;
ceq = [];
end