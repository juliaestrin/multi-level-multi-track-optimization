function result = calculate_Lr(b, w, k_area, N, a, shape, params)

    % Derived geometry
    if strcmp(shape, 'round')
        a         = b;                      % a is ignored for round
        Ac_center = pi * (b/2)^2;
        V_center  = pi*(b/2)^2 * (w/4 + params.h_pcb);
    elseif strcmp(shape, 'rectangle')
        % a passed in as independent variable
        Ac_center = a * b;
        V_center  = a * b * (w/4 + params.h_pcb);
    else
        error('calculate_Lr: unknown shape ''%s''. Use ''round'' or ''rectangle''.', shape);
    end

    Ac_legs     = k_area * Ac_center;
    w_leg       = Ac_legs / a;
    w_tot       = 2*w + b + 2*w_leg;
    w_winding   = w - 2*params.s_ct;
    h_w         = w/4 + params.h_pcb;
    l           = a + 2*w_winding + 2*params.s_ct;
    h           = h_w + 2*w_leg;

    % Core sub-volumes
    V_legs      = 2 * w_leg * a * h_w;
    V_yokes     = 2 * w_tot * a * w_leg;
    Vc          = V_center + V_legs + V_yokes;

    % Initialize result struct
    result.a                = a;
    result.b                = b;
    result.w                = w;
    result.N                = N;
    result.k_area           = k_area;
    result.w_winding        = w_winding;
    result.h_w              = h_w;
    result.w_tot            = w_tot;
    result.l                = l;
    result.h                = h;
    result.Ac_center        = Ac_center;
    result.Ac_legs          = Ac_legs;
    result.w_leg            = w_leg;
    result.Vc               = Vc;
    result.V_center         = V_center;
    result.V_legs           = V_legs;
    result.V_yokes          = V_yokes;
    result.A_footprint      = l * w_tot;
    result.s_ct             = params.s_ct;
    result.h_pcb            = params.h_pcb;
    result.shape            = shape;
    result.gap_loc          = params.gap_loc;

    % Check validity
    if w_winding <= 0 || Ac_center <= 0 || Ac_legs <= 0 || Vc <= 0 || h_w <= 0 || w_leg <= 0
        result = set_invalid(result, 'Negative geometry');
        return;
    end

    % Air gap
    %result.lg = params.u0 * N^2 / params.Lr * (1/(1/2 * 1/Ac_legs + 1/Ac_center));
    result.lg = params.u0 * N^2 / (params.Lr * (1/Ac_center + 1/2*1/Ac_legs)); 
    result.R_legs    = result.lg / (params.u0 * Ac_legs);
    result.R_center  = result.lg / (params.u0 * Ac_center);
    result.R_tot     = 1/2 * result.R_legs + result.R_center;

    % Flux densities
    result.Phi_tot            = N * params.I / result.R_tot;
    result.Bmax_center = result.Phi_tot / Ac_center;
    result.Bmax_legs   = (result.Phi_tot/2) / Ac_legs;

    % Core loss
    % result.P_core = params.k * params.f^(params.alpha) * ( ...
    %     V_center           * result.Bmax_center^(params.beta) + ...
    %     (V_legs + V_yokes) * result.Bmax_legs^(params.beta) );
    % 
    % result.Pv = result.P_core / Vc * 1e-3;

    result.P_core_steinmetz = params.k * params.f^(params.alpha) * ( ...
        V_center * result.Bmax_center^(params.beta) + ...
        (V_legs + V_yokes) * result.Bmax_legs^(params.beta) );

    sigma_1MHz = params.sigma_1MHz; 

    Pv_eddy_center = (8/pi^2) * (pi/6) * sigma_1MHz * b^2 * params.f^2 * result.Bmax_center^2; 
    P_eddy_center = Pv_eddy_center * V_center; 

    Pv_eddy_legs = (8/pi^2) * (pi/6) * sigma_1MHz * w_leg^2 * params.f^2 * result.Bmax_legs^2;
    P_eddy_legs = Pv_eddy_legs * (V_legs + V_yokes); 

    result.P_eddy = P_eddy_center + P_eddy_legs; 

    result.P_core = result.P_core_steinmetz + result.P_eddy; 

    % Winding resistance
    if strcmp(shape, 'round')
        result.Rdc = calculate_Rdc_round( ...
            100, b/2 + params.s_ct, b/2 + params.s_ct + w_winding, ...
            params.sigma_cu, params.t_cu_pri);
    elseif strcmp(shape, 'rectangle')
        result.Rdc = calculate_Rdc_rectangle( ...
            100, b + 2*params.s_ct, b + 2*params.s_ct + w_winding, ...
               a + 2*params.s_ct, a + 2*params.s_ct + w_winding, ...
            params.sigma_cu, params.t_cu_pri);
    end

    % Conduction losses
    result.P_cond = calculate_Pcond_Lr( ...
        params.I, params.f, ...
        result.Rdc, N, ...
        params.t_cu_pri, ...
        params.rho_cu, params.u0);

    result.P_total = result.P_core + result.P_cond;

    if result.P_total < 0 || ~isfinite(result.P_total)
        result = set_invalid(result, 'Invalid losses');
        return;
    end

    result.is_valid    = true;
    result.fail_reason = '';
end

function result = set_invalid(result, reason)
    result.is_valid     = false;
    result.fail_reason  = reason;
    result.Bmax_center  = Inf;
    result.Bmax_legs    = Inf;
    result.lg           = 0;
    result.Rdc          = Inf;
    result.P_core       = Inf;
    result.P_cond       = Inf;
    result.P_total      = Inf;
    result.Pv           = Inf;
    result.Ac_center    = 0;
    result.Ac_legs      = 0;
    result.w_leg        = 0;
    result.V_center     = 0;
    result.V_legs       = 0;
    result.V_yokes      = 0;
    result.Vc           = 0;
    result.h_w          = 0;
    result.w_tot        = 0;
    result.l            = 0;
    result.h            = 0;
end