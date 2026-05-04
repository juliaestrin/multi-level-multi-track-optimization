% OPTIMIZE_VIRT_2LMT Minimize P_total by optimizing a, b, w using fmincon
function [x_opt, P_opt, results] = optimize_VIRT_2LMT(w_max, l_max, params)

    % Bounds [m] - all positive
    lb = [1e-3, 1e-3, 2*params.s_ct];  % [a_min, b_min, w_min]
    ub = [2000e-3, 2000e-3, 2000e-3];  % [a_max, b_max, w_max]
    x0 = [5e-3, 5e-3, 5e-3];          % Initial guess
    
    % Objective: minimize P_total
    objective = @(x) get_Ptotal(x, params);
    
    % Constraints
    nonlcon = @(x) constraints(x, w_max, l_max, params);
    
    % Run fmincon
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off');
    [x_opt, P_opt] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, opts);
    
    % Get full results at optimum
    if strcmp(params.centerpost_shape, 'round')
        results = calculateVIRT_round_2LMT(x_opt(2), x_opt(3), params);
    else
        results = calculateVIRT_square(x_opt(1), x_opt(2), x_opt(3), params);
    end
    
    results.topology = params.topology;
    
    % Compute and store final transformer temperature
    results.T_tx = calculate_transformer_temp( ...
        results.P_total, results.Ac, results.h_w, ...
        params.R_plate, params.Area_plate, params.T_water, ...
        params.sig_grease, params.d_grease);
    end
    
    % -----------------------------------------------------------------------
    
    function P = get_Ptotal(x, params)
    if strcmp(params.centerpost_shape, 'round')
        r = calculateVIRT_round_2LMT(x(2), x(3), params);
    else
        r = calculateVIRT_square(x(1), x(2), x(3), params);
    end
    
    if r.is_valid
        P = r.P_total;
    else
        P = 1e10;  % Penalty for invalid designs
    end
    end
    
    % -----------------------------------------------------------------------
    
    function [c, ceq] = constraints(x, w_max, l_max, params)
    w = x(3);
    s_ct = params.s_ct;
    
    if strcmp(params.centerpost_shape, 'round')
        r = calculateVIRT_round_2LMT(x(2), x(3), params);
    else
        r = calculateVIRT_square(x(1), x(2), x(3), params);
    end
    
    % Geometry
    w_tot     = 2*w + 2*x(2);              % Total width
    w_winding = w - 2*s_ct;
    l         = x(1) + 2*w_winding + 2*s_ct;
    
    % Transformer temperature at this design point
    T_tx = calculate_transformer_temp( ...
        r.P_total, r.Ac, r.h_w, ...
        params.R_plate, params.Area_plate, params.T_water, ...
        params.sig_grease, params.d_grease);
    
    % Inequality constraints: c(i) <= 0
    c(1) = r.Bmax - params.Bsat;   % Bmax <= Bsat
    c(2) = 2*s_ct - w;             % w > 2*s_ct (positive winding width)
    c(3) = -r.Vc;                  % Vc > 0
    c(4) = -r.Ac;                  % Ac > 0
    c(5) = w_tot - w_max;          % w_tot <= w_max
    c(6) = l - l_max;              % l <= l_max
    c(7) = T_tx - params.T_max;    % T_tx <= T_max  <-- NEW
    
    ceq = [];
end