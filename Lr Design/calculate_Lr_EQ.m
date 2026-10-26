function result = calculate_Lr_EQ(b, w, k_area, N, params)
    

    k    = 7.9656e+09;
    beta = 3.3019;

    a           = b;
    Ac_center   = pi * (b/2)^2;
    Ac_legs     = k_area * Ac_center;          % outer leg area (each)
    w_leg       = Ac_legs / a;                 % outer leg width
    w_tot       = 2*w + b + 2*w_leg;
    w_winding   = w - 2*params.s_ct;
    %h_w         = w/4 + params.h_pcb;
    h_w         = 2 * params.h_pcb;
    l           = a + 2*w_winding + 2*params.s_ct;
    h           = h_w + 2*w_leg;

    % Core sub-volumes
    V_center    = pi*(b/2)^2 * h_w;                   % round center post
    V_legs      = 2 * w_leg * a * h_w;                % two outer legs
    V_yokes     = 2 * w_tot * a * w_leg;              % top + bottom yokes
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

    % Air gap
    % result.lg = (N^2 * params.u0 * Ac_center)/params.Lr; 
    % result.R_center = result.lg/(params.u0 * Ac_center); 
    % result.Phi_tot = N*params.I/result.R_center; 
    % result.Bmax_center = result.Phi_tot / Ac_center;
    % result.Bmax_legs   = (result.Phi_tot/2) / Ac_legs;
    %result.h_w = result.h_w + result.lg; % increase the leg height to account for large gap

   
    % result.lg = params.u0 * N^2 / (params.Lr * (1/Ac_center + 1/2*1/Ac_legs)); 
    % result.R_legs    = result.lg / (params.u0 * Ac_legs);
    % result.R_center  = result.lg / (params.u0 * Ac_center);
    % result.R_tot     = 1/2 * result.R_legs + result.R_center;

    % Air gap - solve numerically with fringe effects
    result.lg = solve_lg_fringe(params.Lr, N, params.u0, Ac_center, Ac_legs, b, w_leg, h_w);

% Recompute R_tot using the same fringe-corrected formulas for consistency
lg = result.lg;
denom_c    = params.u0 * (b/(2*lg) + (2/pi)*(1 + log(pi*h_w/(4*lg))));
R_c        = (1/denom_c) * (1/b) * (1/2);

denom_si   = params.u0 * (w_leg/(2*lg) + (2/pi)*(1 + log(pi*h_w/(4*lg))));
R_si       = (1/denom_si) * (1/b);

denom_so   = params.u0 * (w_leg/(2*lg) + (1/pi)*(1 + log(pi*w_leg/(2*lg))));
R_so       = (1/denom_so) * (1/b);

R_s        = (R_si * R_so) / (R_si + R_so);

result.R_center = R_c;
result.R_legs   = R_s;
result.R_tot    = 0.5 * R_s + R_c;


    % 
    % % Flux densities
    result.Phi_tot            = N * params.I / result.R_tot;
    result.Bmax_center = result.Phi_tot / Ac_center;
    result.Bmax_legs   = (result.Phi_tot/2) / Ac_legs;

    % Core loss
    % result.P_core = params.k * params.f^(params.alpha) * ( ...
    %     V_center           * result.Bmax_center^(params.beta) + ...
    %     (V_legs + V_yokes) * result.Bmax_legs^(params.beta) );
    % 
    % result.Pv = result.P_core / Vc * 1e-3;

    % result.P_core_steinmetz = params.k * params.f^(params.alpha) * ( ...
    %     V_center * result.Bmax_center^(params.beta) + ...
    %     (V_legs + V_yokes) * result.Bmax_legs^(params.beta) );


    result.P_core_steinmetz = k * ( ...
        V_center * result.Bmax_center^(beta) + ...
        (V_legs + V_yokes) * result.Bmax_legs^(beta) );

    sigma_1MHz = params.sigma_1MHz; 

   Pv_eddy_center = (8/pi^2) * (pi/4) * sigma_1MHz * Ac_center * params.f^2 * result.Bmax_center^2; 
    P_eddy_center = Pv_eddy_center * V_center; 

    Pv_eddy_legs = (8/pi^2) * (pi/6) * sigma_1MHz * w_leg^2 * params.f^2 * result.Bmax_legs^2;
    P_eddy_legs = Pv_eddy_legs * (V_legs + V_yokes); 
    
    result.P_eddy = P_eddy_center + P_eddy_legs; 

    result.P_core = result.P_core_steinmetz + result.P_eddy; 
    
    result.Pv = result.P_core / Vc * 1e-3;   % [kW/m^3]
    
    result.Rdc = calculate_Rdc_round(100, b/2 + params.s_ct, ...
        b/2 + params.s_ct + w_winding, params.sigma_cu, params.t_cu_pri);
   

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

function lg = solve_lg_fringe(Lr, N, u0, Ac_center, Ac_legs, b, w_leg, h_w)
% Solves for air gap length lg given target inductance Lr,
% using reluctance model with fringe corrections from handwritten equations.
%
% Reluctance model (from notes):
%   R_c         = 1/(u0 * [b/(2*lg) + 2/pi*(1 + ln(pi*h_w/(4*lg)))]) * 1/b * 1/2
%   R_s,inside  = 1/(u0 * [w_leg/(2*lg) + 2/pi*(1 + ln(pi*h_w/(4*lg)))]) * 1/b
%   R_s,outside = 1/(u0 * [w_leg/(2*lg) + 1/pi*(1 + ln(pi*w_leg/(2*lg)))]) * 1/b
%   R_s         = R_s_in * R_s_out / (R_s_in + R_s_out)   [parallel]
%   R_net       = 1/2 * R_s + R_c
%   L           = N^2 / R_net  -->  solve for lg

    function R = R_center(lg)
        % Center post reluctance (round post, gapped on both faces -> *1/2)
        denom = u0 * (b/(2*lg) + (2/pi)*(1 + log(pi*h_w/(4*lg))));
        R = 1/denom * (1/b) * (1/2);
    end

    function R = R_s_inside(lg)
        % Outer leg reluctance seen from winding window (inside face)
        denom = u0 * (w_leg/(2*lg) + (2/pi)*(1 + log(pi*h_w/(4*lg))));
        R = 1/denom * (1/b);
    end

    function R = R_s_outside(lg)
        % Outer leg reluctance seen from outside face
        denom = u0 * (w_leg/(2*lg) + (1/pi)*(1 + log(pi*w_leg/(2*lg))));
        R = 1/denom * (1/b);
    end

    function R = R_legs(lg)
        Ri = R_s_inside(lg);
        Ro = R_s_outside(lg);
        R  = (Ri * Ro) / (Ri + Ro);   % parallel combination
    end

    function R = R_net(lg)
        R = 0.5 * R_legs(lg) + R_center(lg);
    end

% Residual: L(lg) - Lr = 0
% L = N^2 / R_net  -->  N^2/R_net - Lr = 0
residual = @(lg) N^2 / R_net(lg) - Lr;

% Initial bracket: search from 1 um to 10 mm
lg_lo = 1e-6;
lg_hi = 10e-3;

% Safety check: ensure bracket contains a sign change
% if sign(residual(lg_lo)) == sign(residual(lg_hi))
%     warning('solve_lg_fringe: no sign change in [%g, %g]. Falling back to simple model.', lg_lo, lg_hi);
%     % Fallback: original analytical expression (no fringe)
%     lg = u0 * N^2 / (Lr * (1/Ac_center + 0.5/Ac_legs));
%     return;
% end

opts = optimset('TolX', 1e-12, 'Display', 'off');
lg   = fzero(residual, [lg_lo, lg_hi], opts);
end