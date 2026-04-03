function [x_opt, P_opt, results] = optimize_VIRT_v2(w_max, l_max, params)
% OPTIMIZE_VIRT_v2 Minimize P_total by optimizing a, b, w using fmincon

    % Bounds [m] - all positive
    lb = [1e-3, 1e-3, 2*params.s_ct];     % [a_min, b_min, w_min]
    ub = [2000e-3, 2000e-3, 2000e-3];        % [a_max, b_max, w_max]
    x0 = [5e-3, 5e-3, 5e-3];              % Initial guess

    % Objective: minimize P_total
    objective = @(x) get_Ptotal(x, params);

    % Constraints
    nonlcon = @(x) constraints(x, w_max, l_max, params);

    % Run fmincon
    opts = optimoptions('fmincon', 'Algorithm', 'sqp');
    [x_opt, P_opt] = fmincon(objective, x0, [], [], [], [], lb, ub, nonlcon, opts);

    % Get full results at optimum
    
    if strcmp(params.centerpost_shape, 'round')
        results = calculateVIRT_round(x_opt(2), x_opt(3), params);
    else 
        results = calculateVIRT_square(x_opt(1), x_opt(2), x_opt(3), params);
    end 
end

function P = get_Ptotal(x, params)
    
    if strcmp(params.centerpost_shape, 'round')
        r = calculate_VIRT_circularPost_design_v2(x(2), x(3), params);
    else 
        r = calculate_VIRT_design_v2(x(1), x(2), x(3), params);
    end 
    if r.is_valid
        P = r.P_total;
    else
        P = 1e10;  % Penalty for invalid designs
    end
end

function [c, ceq] = constraints(x, w_max, l_max, params)
    a = x(1); 
    b = x(2);
    w = x(3);
    s_ct = params.s_ct;
    
    r = calculate_VIRT_design_v2(x(1), x(2), x(3), params);

    % Inequality constraints: c <= 0
    w_tot = 2*w + 2*b;                  % Total width
    w_winding = w - 2*params.s_ct;
    l = a + 2*w_winding + 2*params.s_ct; 

    c(1) = r.Bmax - params.Bsat;        % Bmax <= Bsat
    c(2) = 2*s_ct - w;                  % w > 2*s_ct (positive winding width)
    c(3) = -r.Vc;                       % Vc > 0 (positive core volume)
    c(4) = -r.Ac;                       % Ac > 0 (positive core area)
    c(5) = w_tot - w_max;              % w_tot <= w_max
    c(6) = l - l_max; 
    ceq = [];
end