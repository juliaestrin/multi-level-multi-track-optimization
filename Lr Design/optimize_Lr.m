function [x_opt, P_opt, results] = optimize_Lr(params)
% optimize_Lr  Minimize total losses of the resonant inductor
%
% Optimization variables: x = [b, w, k_area]
%   b      : center post diameter        [m]
%   w      : winding window width        [m]
%   k_area : ratio Ac_legs / Ac_center   [-]
%   N      : swept as integer externally [-]

N_range = params.N_min : params.N_max;   % integer sweep

% Pre-allocate
P_best    = Inf;
x_best    = [];
N_best    = NaN;

for N = N_range
    [x_N, P_N, flag] = optimize_for_N(N, params);
    if flag > 0 && P_N < P_best      % only accept feasible solutions
        P_best = P_N;
        x_best = x_N;
        N_best = N;
    end
end

% Pack optimal solution
x_opt    = [x_best, N_best];
P_opt    = P_best;
results  = calculate_Lr(x_best(1), x_best(2), x_best(3), N_best, params);

end

% ═════════════════════════════════════════════════════════════════════════
% LOCAL FUNCTIONS
% ═════════════════════════════════════════════════════════════════════════

function [x_opt, P_opt, exitflag] = optimize_for_N(N, params)
% Run fmincon over [b, w, k_area] for a fixed integer N

lb = [1e-3,   2*params.s_ct,  0.1];
ub = [10e-3, 10e-3,         0.5];
x0 = [5e-3,   5e-3,           2];

objective = @(x) get_Ptotal(x, N, params);
nonlcon   = @(x) constraints(x, N, params);

opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');

[x_opt, P_opt, exitflag] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, opts);

end

function P = get_Ptotal(x, N, params)
r = calculate_Lr(x(1), x(2), x(3), N, params);
if r.is_valid
    P = r.P_total;
else
    P = 1e10;
end
end

function [c, ceq] = constraints(x, N, params)
s_ct = params.s_ct;
r    = calculate_Lr(x(1), x(2), x(3), N, params);

% Inequality constraints: c(i) <= 0
c(1) = r.Bmax_center - params.Bsat;   % center post flux density
c(2) = r.Bmax_legs   - params.Bsat;   % outer leg flux density
c(3) = 2*s_ct - x(2);                 % positive winding width (w > 2*s_ct)
c(4) = -r.Vc;                         % Vc > 0
c(5) = -r.Ac_center;                  % Ac_center > 0
c(6) = -r.Ac_legs;                    % Ac_legs > 0

c(~isfinite(c)) = 1e10;               % guard NaN / Inf
ceq = [];
end