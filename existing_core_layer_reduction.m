% Existing Core: Transformer Loss vs. Winding Layer Count
%
% Author: Julia Estrin
%
% Description:
%   Calculates the total transformer loss for the cores we already have
%   (fixed core geometry below) when the PCB winding is changed to reduce
%   the layer count.

clear; close all; clc;

%% ========================= SETUP =========================
addpath('Transformer Design');
addpath('LLC Design');
addpath('Lr Design');

%% ================= EXISTING CORE GEOMETRY =================
b      = 27.7e-3;   % [m]    center post diameter
w      = 8.0e-3;    % [m]    winding window width
k_area = 0.3999;    % [-]    Ac_legs / Ac_center

%% ================= DESIGN PARAMETERS =================
% Values from TwoLevelMutliTrack_design.m (F80 from get_core_material,
% LLC values from designLLC_v3 with the baseline specs)
design_params = struct( ...
    ...   % --- Electrical / LLC ---
    'Vin',      1500,                   ...   % [V]     input voltage
    'f',        1e6,                    ...   % [Hz]    f0
    'I',        13.856324316024601,     ...   % [A]     peak resonant current (Ir_rms = 9.798 A)
    'Lu',       1.9404692383093306e-05, ...   % [H]     magnetizing inductance Lm
    'nt',       2,                      ...   % [-]     number of tracks
    'np',       8,                      ...   % [-]     primary turns per track, N/nt (N = 16)
    ...   % --- Core material (F80) ---
    'k',        6.43007e-05,            ...   % F80 Steinmetz k
    'alpha',    2.15651264669503,       ...   % F80 Steinmetz alpha
    'beta',     2.549625682,            ...   % F80 Steinmetz beta
    'sigma_1MHz', 1.363274,             ...   % [S/m]   F80 conductivity
    ...   % --- Core shape and gap ---
    'centerpost_shape', 'round',        ...
    'gap_loc',  'all',                  ...   % gap in all legs
    'u0',       4*pi*1e-7,              ...   % [H/m]   permeability of free space
    ...   % --- Winding / PCB ---
    'stackup',  '5Layer',               ...   % unused label in calculate_Pcond_2LMT
    'h_pcb',    0.002108,               ...   % [m]     thickness of pcb
    's_ct',     1e-3,                   ...   % [m]     copper to core distance
    's_tt_pri', 0.1238e-3,              ...   % [m]     primary turn-to-turn spacing (93.75 V turn-to-turn)
    't_cu_pri', 2*35e-6,                ...   % [m]     2 oz
    't_cu_sec', 2*35e-6,                ...   % [m]     2 oz
    'rho_cu',   2.2e-8,                 ...   % [Ohm*m] copper resistivity
    'sigma_cu', 1/2.2e-8,               ...   % [S/m]   copper conductivity
    ...   % --- Thermal (cold plate + grease) ---
    'R_plate',  0.08,                   ...   % [C/W]
    'Area_plate', 152.4e-3 * 76.2e-3,   ...   % [m^2]
    'T_water',  45,                     ...   % [C]
    'sig_grease', 4,                    ...   % [W/(mK)]
    'd_grease', 0.127e-3                ...   % [m]
    );

TX_design = calculateVIRT_6layer_S_4Px2T_S(b, w, k_area, design_params);

% Transformer temperature (same model as optimize_VIRT_2LMT)
TX_design.T_tx = calculate_transformer_temp( ...
    TX_design.P_total, TX_design.Ac_center, TX_design.h_w, ...
    design_params.R_plate, design_params.Area_plate, design_params.T_water, ...
    design_params.sig_grease, design_params.d_grease);

material = get_core_material('F80');
VIRT3Dfigure_round(TX_design, material, TX_design.T_tx, 0);

function result = calculateVIRT_6layer_S_4Px2T_S(b, w, k_area, params)
    result = struct();

    % Derived geometry
    a           = b;
    Ac_center   = pi * (b/2)^2;
    Ac_legs     = k_area * Ac_center;          % outer leg area (each)
    w_leg       = Ac_legs / a;                 % outer leg width
    w_tot       = 2*w + b + 2*w_leg;
    w_winding   = w - 2*params.s_ct;
    h_w         = w/4 + params.h_pcb;
    l           = a + 2*w_winding + 2*params.s_ct;
    h           = h_w + 2*w_leg;
    
    % Core sub-volumes
    V_center    = pi*(b/2)^2 * h_w;                   % round center post
    V_legs      = 2 * w_leg * a * h_w;                % two outer legs
    V_yokes     = 2 * w_tot * a * w_leg;              % top + bottom yokes
    Vc          = V_center + V_legs + V_yokes;
    
    
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
    w_turn_pri  = (w_winding - params.s_tt_pri) / 2;     % radial width of each primary turn

    % Primary turn radii (inner turn, gap s_tt_pri, outer turn)
    r_in_pri_inner  = b/2 + params.s_ct;                  % inner turn, inner radius
    r_out_pri_inner = r_in_pri_inner + w_turn_pri;        % inner turn, outer radius
    r_in_pri_outer  = r_out_pri_inner + params.s_tt_pri;  % outer turn, inner radius
    r_out_pri_outer = r_in_pri_outer + w_turn_pri;        % outer turn, outer radius (= b/2 + s_ct + w_winding)

    result.Rdc_pri_inner = calculate_Rdc_round(r_in_pri_inner, r_out_pri_inner, params.sigma_cu, params.t_cu_pri);
    result.Rdc_pri_outer = calculate_Rdc_round(r_in_pri_outer, r_out_pri_outer, params.sigma_cu, params.t_cu_pri);

    % Secondary turn radii (one full-width turn)
    r_in_sec  = b/2 + params.s_ct;                        % inner radius
    r_out_sec = r_in_sec + w_winding;                     % outer radius

    result.Rdc_sec = calculate_Rdc_round(r_in_sec, r_out_sec, params.sigma_cu, params.t_cu_sec);
    
    % Conduction losses
    [result.P_pri, result.P_sec] = calculate_Pcond_6layer_S_4Px2T_S( ...
        params.I, params.f, ...
        result.Rdc_pri_inner, result.Rdc_pri_outer,result.Rdc_sec, ...
        params.t_cu_pri, params.t_cu_sec, ...
        params.rho_cu, params.u0);
    
    result.P_total = result.P_core + result.P_pri + result.P_sec;
    
    %% Initialize result struct with geometry
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

    %% Check for valid loss values
    if result.P_total < 0 || ~isfinite(result.P_total)
        result = set_invalid(result, 'Invalid losses');
        return;
    end
    
    % Valid design
    result.is_valid    = true;
    result.fail_reason = '';
end

function [P_pri, P_sec] = calculate_Pcond_6layer_S_4Px2T_S(I, f, Rdc_pri_inner, Rdc_pri_outer, Rdc_sec, t_cu_pri, t_cu_sec, rho_cu, u0)

    %% Skin Depth and Normalized Copper Thickness
    skin_depth  = sqrt(rho_cu / (pi * f * u0));
    delta_pri   = t_cu_pri / skin_depth;
    delta_sec   = t_cu_sec / skin_depth;

    %% Dowell Coefficients
    F_pri = (sinh(2*delta_pri) + sin(2*delta_pri)) / (cosh(2*delta_pri) - cos(2*delta_pri));
    G_pri = (sinh(delta_pri)   - sin(delta_pri))   / (cosh(delta_pri)   + cos(delta_pri));
    
    F_sec = (sinh(2*delta_sec) + sin(2*delta_sec)) / (cosh(2*delta_sec) - cos(2*delta_sec));
    G_sec = (sinh(delta_sec)   - sin(delta_sec))   / (cosh(delta_sec)   + cos(delta_sec));


    %% MMF at each layer boundary [x I]
    %        top   S1   P1   P2   P3   P4   S2   bottom
    MMF = [   0,   -4,  -2,   0,   2,   4,   0];
    N_pri = 2;                                   % primary turns per layer

    %% Secondary layers (1 full-width turn each)
    P_S1 = 1/2 * Rdc_sec * delta_sec * ((MMF(1) - MMF(2))^2 * F_sec + 2 * MMF(1) * MMF(2) * G_sec) * I^2;
    P_S2 = 1/2 * Rdc_sec * delta_sec * ((MMF(6) - MMF(7))^2 * F_sec + 2 * MMF(6) * MMF(7) * G_sec) * I^2;

    %% Primary layers (inner and outer turn of each layer)
    P_P1_inner = 1/2 * Rdc_pri_inner / N_pri^2 * delta_pri * ((MMF(2) - MMF(3))^2 * F_pri + 2 * MMF(2) * MMF(3) * G_pri) * I^2;
    P_P1_outer = 1/2 * Rdc_pri_outer / N_pri^2 * delta_pri * ((MMF(2) - MMF(3))^2 * F_pri + 2 * MMF(2) * MMF(3) * G_pri) * I^2;

    P_P2_inner = 1/2 * Rdc_pri_inner / N_pri^2 * delta_pri * ((MMF(3) - MMF(4))^2 * F_pri + 2 * MMF(3) * MMF(4) * G_pri) * I^2;
    P_P2_outer = 1/2 * Rdc_pri_outer / N_pri^2 * delta_pri * ((MMF(3) - MMF(4))^2 * F_pri + 2 * MMF(3) * MMF(4) * G_pri) * I^2;

    P_P3_inner = 1/2 * Rdc_pri_inner / N_pri^2 * delta_pri * ((MMF(4) - MMF(5))^2 * F_pri + 2 * MMF(4) * MMF(5) * G_pri) * I^2;
    P_P3_outer = 1/2 * Rdc_pri_outer / N_pri^2 * delta_pri * ((MMF(4) - MMF(5))^2 * F_pri + 2 * MMF(4) * MMF(5) * G_pri) * I^2;

    P_P4_inner = 1/2 * Rdc_pri_inner / N_pri^2 * delta_pri * ((MMF(5) - MMF(6))^2 * F_pri + 2 * MMF(5) * MMF(6) * G_pri) * I^2;
    P_P4_outer = 1/2 * Rdc_pri_outer / N_pri^2 * delta_pri * ((MMF(5) - MMF(6))^2 * F_pri + 2 * MMF(5) * MMF(6) * G_pri) * I^2;

    %% Totals
    P_sec = P_S1 + P_S2;
    P_pri = P_P1_inner + P_P1_outer + P_P2_inner + P_P2_outer + ...
            P_P3_inner + P_P3_outer + P_P4_inner + P_P4_outer;

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