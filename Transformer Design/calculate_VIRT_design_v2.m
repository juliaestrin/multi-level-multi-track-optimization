% Generate an EQ Core Design for a 1/2-Turn VIRT
%
% Author:  Julia Estrin, Qijia Li
% Date:    03-18-2026
%
% Description:
%   Designs an EQ-core planar transformer for a 1/2-turn VIRT topology.
%   Given a maximum core loss density and window height, this function
%   computes all core and winding geometry, then evaluates core loss and
%   AC copper loss (via Dowell's method) to return a fully characterized
%   design result.
%
%   Core geometry assumptions:
%     - Cylindrical center post (EQ core style)
%     - Outer leg area = 1/2 center post area  =>  uniform B across all sections
%     - Yoke area = leg area  =>  h_yoke = l_leg
%     - Square winding footprint (l_winding == w_winding) for this topology
%     - Air gaps in outer legs (two gaps in series)
%
%   Magnetic circuit model (updated 02-27-2026):
%     - Square wave voltage Vso = Np * nt * Vo generated across primaries
%     - Peak triangular flux: Ac = Vso / (4 * f * Np * Bmax)
%     - Air gap per leg: lg = u0 * Np^2 * Ac / (nt * Lu)
%     - Factor nt accounts for tracks in series
%
% Inputs:
%   Pv_max    - Maximum allowable core loss density [W/m^3]
%   w_height  - Winding window height [m]
%   params    - Parameter struct with fields:
%                 .k        - Steinmetz coefficient k  (P = k * B^beta)
%                 .beta     - Steinmetz exponent beta
%                 .Lu       - Inductance requirement per track [H]
%                 .I        - Peak/RMS winding current [A]
%                 .Np       - Number of primary turns
%                 .nt       - Number of secondary tracks
%                 .Vo       - Output voltage per track [V]
%                 .w_b      - Winding breadth (window width) [m]
%                 .f        - Switching frequency [Hz]
%                 .sigma_cu - Copper electrical conductivity [S/m]
%                 .rho_cu   - Copper electrical resistivity [Ohm·m]
%                 .t_cu_pri - Primary copper trace thickness [m]
%                 .t_cu_sec - Secondary copper trace thickness [m]
%                 .u0       - Permeability of free space [H/m]
%                 .stackup  - Winding stackup configuration string
%
% Outputs:
%   result    - Struct containing all design variables and losses:
%                 .Bmax         - Peak flux density [T]
%                 .Ac           - Required center leg cross-section area [m^2]
%                 .r_centerpost - Center post radius [m]
%                 .w_core       - Core width (= center post diameter) [m]
%                 .l_leg        - Outer leg length [m]
%                 .h_yoke       - Yoke height (= l_leg) [m]
%                 .l_core       - Total core length [m]
%                 .h_core       - Total core height [m]
%                 .w_b          - Winding window breadth/depth [m]
%                 .l_winding    - Winding outer dimension, length [m]
%                 .w_winding    - Winding outer dimension, width [m]
%                 .V_total      - Total core volume [m^3]
%                 .A_footprint  - PCB footprint area [m^2]
%                 .lg           - Air gap length per outer leg [m]
%                 .P_core       - Core loss [W]
%                 .Rdc_pri      - Primary DC winding resistance [Ohm]
%                 .Rdc_sec      - Secondary DC winding resistance [Ohm]
%                 .P_pri        - Primary AC conduction loss [W]
%                 .P_sec        - Secondary AC conduction loss [W]
%                 .P_total      - Total loss (core + primary + secondary) [W]
%                 .Cps          - Primary-to-secondary parasitic capacitance [F]
%
% Usage Example:
%   params.k        = 1.5e-3;
%   params.beta     = 2.2;
%   params.Lu       = 10e-6;
%   params.I        = 10;
%   params.Np       = 4;
%   params.nt       = 2;
%   params.Vo       = 48;
%   params.w_b      = 5e-3;
%   params.f        = 100e3;
%   params.sigma_cu = 5.8e7;
%   params.rho_cu   = 1.72e-8;
%   params.t_cu_pri = 35e-6;
%   params.t_cu_sec = 35e-6;
%   params.u0       = 4*pi*1e-7;
%   params.stackup  = '3layer';
%   result = calculate_VIRT_design(5e4, 3e-3, params);

function result = calculate_VIRT_design_v2(Pv_max, w_height, h_core_max, params)

    %% Flux Density
    Bmax = (Pv_max / (params.k * params.f ^ params.alpha))^(1 / params.beta);

    %% Core Geometry
    Vso = params.Np * params.nt * params.Vo;

    % Required center leg area
    Ac = Vso / (4 * params.f * params.Np * Bmax);

    % Air gap
    lg = params.u0 * params.Np^2 * Ac / (params.nt * params.Lu);

    % Square centerpost
    w_core = sqrt(Ac);

    % Outer legs
    l_leg  = 0.5 * Ac / w_core;

    % Yoke
    h_yoke = l_leg;

    % Core dimensions
    l_core = 2 * w_height + 2 * l_leg + w_core;
    h_core = 2 * h_yoke + params.w_b;

    %% Absolute height constraint 
    if h_core > h_core_max
        result = struct( ...
            'is_valid',      false, ...
            'fail_reason',   'h_core exceeds h_core_max', ...
            'Bmax',          NaN, ...
            'Ac',            NaN, ...
            'w_core',        NaN, ...
            'l_leg',         NaN, ...
            'h_yoke',        NaN, ...
            'l_core',        NaN, ...
            'h_core',        h_core, ...
            'w_b',           params.w_b, ...
            'l_winding',     NaN, ...
            'w_winding',     NaN, ...
            'V_total',       NaN, ...
            'A_footprint',   NaN, ...
            'lg',            NaN, ...
            'P_core',        NaN, ...
            'Rdc_pri',       NaN, ...
            'Rdc_sec',       NaN, ...
            'P_pri',         NaN, ...
            'P_sec',         NaN, ...
            'P_total',       NaN ...
        );
        return;
    end

    %% Winding Geometry
    l_winding = w_core + 2 * w_height;
    w_winding = w_core + params.w_scale * 2 * w_height;

    %% Core Areas
    A_centerpost = w_core^2;
    A_leg        = l_leg  * w_core;
    A_yoke       = h_yoke * w_core;
    A_footprint  = w_winding * l_core;

    %% Core Volumes
    V_centerpost = A_centerpost * params.w_b;
    V_legs       = 2 * A_leg  * params.w_b;
    V_yokes      = 2 * A_yoke * l_core;
    V_total      = V_centerpost + V_legs + V_yokes;

    %% Core Loss
    P_core = V_total * params.k * params.f^params.alpha * Bmax^params.beta;

    %% DC Resistance
    Rdc_pri = calculate_Rdc_summation(100, ...
        w_core, l_winding, ...
        w_core, w_winding, ...
        params.sigma_cu, params.t_cu_pri);

    Rdc_sec = calculate_Rdc_summation(100, ...
        w_core, l_winding, ...
        w_core, w_winding, ...
        params.sigma_cu, params.t_cu_sec);

    %% AC Copper Loss
    [P_pri, P_sec] = calculate_Pcond( ...
        params.I, params.f, ...
        Rdc_pri, Rdc_sec, ...
        params.t_cu_pri, params.t_cu_sec, ...
        params.rho_cu, params.u0, params.stackup);

    %% Total Loss
    P_total = P_core + P_pri + P_sec;

    %% Package Results
    result = struct( ...
        'is_valid',      true, ...
        'fail_reason',   '', ...
        'Bmax',          Bmax, ...
        'Ac',            Ac, ...
        'w_core',        w_core, ...
        'l_leg',         l_leg, ...
        'h_yoke',        h_yoke, ...
        'l_core',        l_core, ...
        'h_core',        h_core, ...
        'w_b',           params.w_b, ...
        'l_winding',     l_winding, ...
        'w_winding',     w_winding, ...
        'V_total',       V_total, ...
        'A_footprint',   A_footprint, ...
        'lg',            lg, ...
        'P_core',        P_core, ...
        'Rdc_pri',       Rdc_pri, ...
        'Rdc_sec',       Rdc_sec, ...
        'P_pri',         P_pri, ...
        'P_sec',         P_sec, ...
        'P_total',       P_total ...
    );

end