function result = calculateVIRT_roundEQ_2LMT(b, w, k_area, params)
    % Derived geometry
    a           = b;
    Ac_center   = pi * (b/2)^2;
    Ac_legs     = k_areSa * Ac_center;          % outer leg area (each)
    %w_leg       = w + sqrt(2*Ac_legs/pi + w^2);                 % outer leg width
    % t = acos((b/2)/(b/2+w)); 
    % w_leg = (-b * t - 2* t*w + sqrt(t*(4*Ac_legs + t * (b + 2* w)^2)))/(2* t); 

    w_leg = Ac_legs / a; 
    w_tot       = 2*w + b + 2*w_leg;
    w_winding   = w - 2*params.s_ct;
    h_w         = w/4 + params.h_pcb;
    l           = a + 2*params.s_ct + 2*w_winding;
    h           = h_w + 2*w_leg;
    
    % Core sub-volumes
    V_center    = pi*(b/2)^2 * h_w;                   % round center post
    V_legs      = 2 * Ac_legs * h_w;                % two outer legs
    V_yokes     = 2* w_leg * a * w_tot ; % this is approximate             % top + bottom yokes
    Vc          = V_center + V_legs + V_yokes;
    
    % Initialize result struct with geometry
    result.a                = a;
    result.b                = b;
    result.w                = w;
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
    result.centerpost_shape = params.centerpost_shape;
    result.gap_loc          = params.gap_loc;
    
    % Check validity
    if w_winding <= 0 || Ac_center <= 0 || Ac_legs <= 0 || Vc <= 0 || h_w <= 0 || w_leg <= 0
        result = set_invalid(result, 'Negative geometry');
        return;
    end
    
    % Air gap length — depends on gap location
    switch params.gap_loc
        case 'center'
            % Gap only in center post
            result.lg = params.u0 * params.np^2 * Ac_center / (params.nt * params.Lu);
    
        case 'all'
            % Equal gap in all legs: lg_center = lg_leg = lg
            result.lg = params.u0 * params.np^2 / (params.nt * params.Lu * (0.5*1/Ac_legs + 1/Ac_center));
    
        otherwise  
            error('calculateVIRT_round_2LMT: unknown gap_loc ''%s''. Use ''center'' or ''all''.', params.gap_loc);
    end
    
    % Flux densities — center post carries full flux, each outer leg carries half
    % result.Bmax_center = 2*params.Vo / (4*params.f * Ac_center);
    % result.Bmax_legs   = 2*params.Vo / (4*params.f * 2 * Ac_legs);   % legs + yokes
    
    result.Bmax_center = params.Vin / (8 * 8 *params.f * Ac_center);
    result.Bmax_legs   = params.Vin / (8 * 8 *params.f * 2 * Ac_legs);   % legs + yokes

    % Steinmetz Core loss — each region at its own Bmax
    result.P_core_steinmetz = params.k * params.f^(params.alpha) * ( ...
        V_center * result.Bmax_center^(params.beta) + ...
        (V_legs + V_yokes) * result.Bmax_legs^(params.beta) );
    
    sigma_1MHz = params.sigma_1MHz; 

    Pv_eddy_center = (8/pi^2) * (pi/4) * sigma_1MHz * Ac_center * params.f^2 * result.Bmax_center^2; 
    P_eddy_center = Pv_eddy_center * V_center; 

    Pv_eddy_legs = (8/pi^2) * (pi/6) * sigma_1MHz * w_leg^2 * params.f^2 * result.Bmax_legs^2;
    P_eddy_legs = Pv_eddy_legs * (V_legs + V_yokes); 
    
    result.P_eddy = P_eddy_center + P_eddy_legs; 

    result.P_core = result.P_core_steinmetz + result.P_eddy; 
    
    result.Pv = result.P_core / Vc * 1e-3;   % [kW/m^3]
    
    % Winding resistance
    result.Rdc_pri = calculate_Rdc_round(100, b/2 + params.s_ct, ...
        b/2 + params.s_ct + w_winding, params.sigma_cu, params.t_cu_pri);
    result.Rdc_sec = calculate_Rdc_round(100, b/2 + params.s_ct, ...
        b/2 + params.s_ct + w_winding, params.sigma_cu, params.t_cu_sec);
    
    % Conduction losses
    [result.P_pri, result.P_sec] = calculate_Pcond_2LMT( ...
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
    result.is_valid    = true;
    result.fail_reason = '';
    end
    
    function result = set_invalid(result, reason)
    result.is_valid     = false;
    result.fail_reason  = reason;
    result.Bmax_center  = Inf;
    result.Bmax_legs    = Inf;
    result.lg           = 0;
    result.Rdc_pri      = Inf;
    result.Rdc_sec      = Inf;
    result.P_core       = Inf;
    result.P_pri        = Inf;
    result.P_sec        = Inf;
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