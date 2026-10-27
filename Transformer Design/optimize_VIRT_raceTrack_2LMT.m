% function [x_opt, P_opt, result] = optimize_VIRT_raceTrack_2LMT(w_max, l_max, params, x0)
% 
% if nargin < 4
%     x0 = [5e-3, 10e-3, 5e-3];   % [b, a, w]
% end
% 
% lb = [1e-3, 1e-3, 2*params.s_ct];
% ub = [2000e-3, 2000e-3, 2000e-3];
% 
% objective = @(x) get_Ptotal_raceTrack(x, params);
% nonlcon   = @(x) constraints_raceTrack(x, w_max, l_max, params);
% 
% opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
% 
% [x_opt, P_opt] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, opts);
% 
% if isinf(P_opt) || P_opt >= 1e10
%     warning('optimize_VIRT_raceTrack_2LMT: no valid solution found for k_area = %.3f', params.k_area);
% end
% 
% b_opt = x_opt(1);
% a_opt = x_opt(2);
% w_opt = x_opt(3);
% 
% result          = calculateVIRT_raceTrack_2LMT(a_opt, b_opt, w_opt, params.k_area, params);
% result.topology = params.topology;
% result.k_area   = params.k_area;
% result.T_tx     = calculate_transformer_temp( ...
%     result.P_total, result.Ac_center, result.h_w, ...
%     params.R_plate, params.Area_plate, params.T_water, ...
%     params.sig_grease, params.d_grease);
% 
% end
% 
% % -------------------------------------------------------------------------
% function P = get_Ptotal_raceTrack(x, params)
% 
% b = x(1);
% a = x(2);
% w = x(3);
% 
% r = calculateVIRT_raceTrack_2LMT(a, b, w, params.k_area, params);
% 
% if r.is_valid
%     P = r.P_total;
% else
%     P = 1e10;
% end
% 
% end
% 
% % -------------------------------------------------------------------------
% function [c, ceq] = constraints_raceTrack(x, w_max, l_max, params)
% 
% b = x(1);
% a = x(2);
% w = x(3);
% 
% r = calculateVIRT_raceTrack_2LMT(a, b, w, params.k_area, params);
% 
% T_tx = calculate_transformer_temp( ...
%     r.P_total, r.Ac_center, r.h_w, ...
%     params.R_plate, params.Area_plate, params.T_water, ...
%     params.sig_grease, params.d_grease);
% 
% c(1) = r.Bmax_center - params.Bsat;   % flux density limit, center post
% c(2) = r.Bmax_legs   - params.Bsat;   % flux density limit, outer legs
% c(3) = 2*params.s_ct - w;             % minimum winding thickness
% c(4) = -r.Vc;                         % core volume must be positive
% c(5) = -r.Ac_center;                  % center post area must be positive
% c(6) = -r.Ac_legs;                    % leg area must be positive
% c(7) = r.w_tot       - w_max;         % total width fits in footprint
% c(8) = r.l           - l_max;         % total length fits in footprint
% c(9) = T_tx          - params.T_max;  % temperature limit
% c(10) = b - a;                        % racetrack validity: a >= b
% 
% ceq = [];
% 
% end

function [x_opt, P_opt, result] = optimize_VIRT_raceTrack_2LMT(w_max, l_max, params, x0)

% Pre-compute a conservative w_leg estimate for linear constraint tightening.
% w_leg = k_area * Ac_center / a, and Ac_center grows with b and a,
% so we use a small representative geometry to get a lower-bound on w_leg.
% Instead, we simply reduce the effective w_max available to [b, w] by
% reserving space for 2*w_leg using a worst-case estimate.
% The nonlinear constraint is the true enforcer; we just need a reasonable
% starting region.
w_leg_est = params.k_area * (pi*(w_max/8)^2) / (l_max/4) ;
w_max_eff = w_max - 2*w_leg_est;
w_max_eff = max(w_max_eff, 5e-3);   % floor to keep problem feasible

if nargin < 4
    x0 = [w_max_eff/4, l_max/4, w_max_eff/6];   % [b, a, w]
end

% --- Bounds ---
lb = [1e-3,       1e-3,  2*params.s_ct];
ub = [w_max_eff,  l_max, w_max_eff/2  ];

% --- Linear constraints ---
% b + 2*w <= w_max_eff  (conservative: leaves room for 2*w_leg)
% a + 2*w <= l_max + 2*s_ct  (exact for l)
A_lin = [ 1  0  2;
          0  1  2];
b_lin = [ w_max_eff;
          l_max + 2*params.s_ct];

objective = @(x) get_Ptotal_raceTrack(x, params);
nonlcon   = @(x) constraints_raceTrack(x, w_max, l_max, params);

opts = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...
    'Display',   'off');

[x_opt, P_opt] = fmincon(objective, x0, A_lin, b_lin, [], [], lb, ub, nonlcon, opts);

if isinf(P_opt) || P_opt >= 1e10
    warning('optimize_VIRT_raceTrack_2LMT: no valid solution found for k_area = %.3f', params.k_area);
end

b_opt = x_opt(1);
a_opt = x_opt(2);
w_opt = x_opt(3);

result          = calculateVIRT_raceTrack_2LMT(a_opt, b_opt, w_opt, params.k_area, params);
result.topology = params.topology;
result.k_area   = params.k_area;
result.T_tx     = calculate_transformer_temp( ...
    result.P_total, result.Ac_center, result.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);

end

% -------------------------------------------------------------------------
function P = get_Ptotal_raceTrack(x, params)

b = x(1);
a = x(2);
w = x(3);

r = calculateVIRT_raceTrack_2LMT(a, b, w, params.k_area, params);

if r.is_valid
    P = r.P_total;
else
    P = 1e10;
end

end

% -------------------------------------------------------------------------
function [c, ceq] = constraints_raceTrack(x, w_max, l_max, params)

b = x(1);
a = x(2);
w = x(3);

r = calculateVIRT_raceTrack_2LMT(a, b, w, params.k_area, params);

if ~r.is_valid
    c   = 1e6 * ones(10, 1);
    ceq = [];
    return;
end

T_tx = calculate_transformer_temp( ...
    r.P_total, r.Ac_center, r.h_w, ...
    params.R_plate, params.Area_plate, params.T_water, ...
    params.sig_grease, params.d_grease);

c(1)  = r.Bmax_center - params.Bsat;
c(2)  = r.Bmax_legs   - params.Bsat;
c(3)  = 2*params.s_ct - w;
c(4)  = -r.Vc;
c(5)  = -r.Ac_center;
c(6)  = -r.Ac_legs;
c(7)  = r.w_tot       - w_max;     % true enforcer — uses actual w_leg
c(8)  = r.l           - l_max;
c(9)  = T_tx          - params.T_max;
c(10) = b             - a;

ceq = [];

end