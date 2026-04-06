function result = calculateSPVIRT_square(a, b, w, params)
    % Derived geometry
    w_winding = w - 2*params.s_ct;
    h_w = w/4 + params.h_pcb;
    w_tot = 4*w + 3*b;
    l = a + 2*w_winding + 2*params.s_ct;
    h = h_w + b;
    Ac = a * b;
    Vc = w_tot * a * h - 4 * a * w * h_w;

    % Initialize result struct with geometry (always included)
    result.a = a;
    result.b = b;
    result.w = w;
    result.w_winding = w_winding;
    result.h_w = h_w;
    result.w_tot = w_tot;
    result.l = l;
    result.h = h;
    result.Ac = Ac;
    result.Vc = Vc;
    result.A_footprint = l * w_tot;
    result.s_ct = params.s_ct;
    result.h_pcb = params.h_pcb;

    % Check validity - all quantities must be positive
    if w_winding <= 0 || Ac <= 0 || Vc <= 0 || h_w <= 0
        result = set_invalid(result, 'Negative geometry');
        return;
    end

    % Calculate electromagnetic quantities
    result.lg = params.u0 * params.np^2 * Ac / (8 * params.Lu);
    result.Bmax = 2*params.Vo / (4*params.f*Ac);
    result.P_core = Vc * params.k * params.f^(params.alpha) * result.Bmax^(params.beta);
    result.Pv = result.P_core / Vc * 1e-3; % [kW/m^3]

    % Winding resistance for a single phase
    result.Rdc_pri = calculate_Rdc_rectangle(100, ...
            b + 2*params.s_ct, b + 2*params.s_ct + 2*w_winding, ...
            a + 2*params.s_ct, l, ...
            params.sigma_cu, params.t_cu_pri);

    result.Rdc_sec = calculate_Rdc_rectangle(100, ...
            b + 2*params.s_ct, b + 2*params.s_ct + 2*w_winding, ...
            a + 2*params.s_ct, l, ...
            params.sigma_cu, params.t_cu_sec);

    % Conduction losses
    [result.P_pri, result.P_sec] = calculate_Pcond_SPVIRT( ...
        params.I, params.f, ...
        result.Rdc_pri, result.Rdc_sec, ...
        params.t_cu_pri, params.t_cu_sec, ...
        params.rho_cu, params.u0, params.stackup);

    result.P_total = result.P_core + result.P_pri + result.P_sec;

    % Check for valid loss values
    if result.P_total < 0 || ~isfinite(result.P_total)
        result = set_invalid(result, 'Invalid losses');
        return;
    end

    % Valid design
    result.is_valid = true;
    result.fail_reason = '';
end

function result = set_invalid(result, reason)
    result.is_valid = false;
    result.fail_reason = reason;
    result.Bmax = Inf;
    result.lg = 0;
    result.Rdc_pri = Inf;
    result.Rdc_sec = Inf;
    result.P_core = Inf;
    result.P_pri = Inf;
    result.P_sec = Inf;
    result.P_total = Inf;
    result.Pv = Inf;
end