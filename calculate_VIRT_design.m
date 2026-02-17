% Generate an EQ Core Design for a 1/2-Turn VIRT
%
% Author:  Julia Estrin
% Date:    02-12-2026
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
%
% Inputs:
%   Pv_max    - Maximum allowable core loss density [W/m^3]
%   w_height  - Winding window height [m]
%   params    - Parameter struct with fields:
%                 .k        - Steinmetz coefficient k  (P = k * B^beta)
%                 .beta     - Steinmetz exponent beta
%                 .Lu       - Inductance requirement [H]
%                 .I        - Peak/RMS winding current [A]
%                 .Np       - Number of primary turns
%                 .w_b      - Winding breadth (window width) [m]
%                 .f        - Switching frequency [Hz]
%                 .sigma_cu - Copper electrical conductivity [S/m]
%                 .rho_cu   - Copper electrical resistivity [Ohm·m]
%                 .t_cu_pri - Primary copper trace thickness [m]
%                 .t_cu_sec - Secondary copper trace thickness [m]
%                 .u0       - Permeability of free space [H/m]
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
%                 .P_core       - Core loss [W]
%                 .Rdc_pri      - Primary DC winding resistance [Ohm]
%                 .Rdc_sec      - Secondary DC winding resistance [Ohm]
%                 .P_pri        - Primary AC conduction loss [W]
%                 .P_sec        - Secondary AC conduction loss [W]
%                 .P_total      - Total loss (core + primary + secondary) [W]
%
% Usage Example:
%   params.k        = 1.5e-3;
%   params.beta     = 2.2;
%   params.Lu       = 10e-6;
%   params.I        = 10;
%   params.Np       = 4;
%   params.w_b      = 5e-3;
%   params.f        = 100e3;
%   params.sigma_cu = 5.8e7;
%   params.rho_cu   = 1.72e-8;
%   params.t_cu_pri = 35e-6;
%   params.t_cu_sec = 35e-6;
%   params.u0       = 4*pi*1e-7;
%   result = calculate_VIRT_design(5e4, 3e-3, params);

function result = calculate_VIRT_design(Pv_max, w_height, params)

    %% Flux Density
    Bmax = (Pv_max / params.k)^(1 / params.beta);

    %% Core Geometry
    % Center post sized to satisfy inductance requirement at Bmax
    Ac            = params.Lu * params.I / (params.Np * Bmax);     % [m^2] required center leg area
    r_centerpost  = sqrt(Ac / pi);                                  % [m]   center post radius (cylindrical)
    w_core        = 2 * r_centerpost;                               % [m]   core width = center post diameter

    % Outer legs: area = 1/2 center post area => uniform B in all sections
    l_leg         = 0.5 * Ac / w_core;                             % [m]   outer leg length

    % Yoke: area matches leg area => h_yoke = l_leg
    h_yoke        = l_leg;                                          % [m]   yoke height

    % Overall core envelope
    l_core        = 2 * w_height + 2 * l_leg + 2 * r_centerpost;  % [m]   total core length
    h_core        = 2 * h_yoke + params.w_b;                       % [m]   total core height

    %% Winding Geometry
    % Square winding footprint for this topology (l_winding == w_winding)
    l_winding     = 2 * r_centerpost + 2 * w_height;               % [m]   winding outer length
    w_winding     = 2 * r_centerpost + 1/2 * 2 * w_height;               % [m]   winding outer width

    %% Core Cross-Section Areas
    A_centerpost  = pi * r_centerpost^2;                            % [m^2]
    A_leg         = l_leg  * w_core;                                % [m^2]
    A_yoke        = h_yoke * w_core;                                % [m^2]
    A_footprint   = w_core * l_core;                                % [m^2] PCB footprint

    %% Core Section Volumes
    V_centerpost  = A_centerpost * params.w_b;                      % [m^3] cylindrical center post
    V_legs        = 2 * A_leg    * params.w_b;                      % [m^3] left and right outer legs
    V_yokes       = 2 * A_yoke   * l_core;                          % [m^3] top and bottom yokes
    V_total       = V_centerpost + V_legs + V_yokes;                % [m^3] total core volume

    %% Core Loss
    % B field is uniform across all sections (leg/yoke area = 1/2 center post area)
    P_core        = V_total * params.k * Bmax^params.beta;          % [W]

    %% DC Winding Resistance
    % calculate_Rdc_summation(n, a_0, a_n, b_0, b_n, sigma_cu, t_cu)
    Rdc_pri = calculate_Rdc_summation(100, ...
        w_core, l_winding, ...
        w_core, w_winding, ...
        params.sigma_cu, params.t_cu_pri);

    Rdc_sec = calculate_Rdc_summation(100, ...
        w_core, l_winding, ...
        w_core, w_winding, ...
        params.sigma_cu, params.t_cu_sec);

    %% AC Copper Loss (Dowell's Method)
    % calculate_Pcond_3Layer(I, f, Rdc_pri, Rdc_sec, t_cu_pri, t_cu_sec, rho_cu, u0)
    [P_pri, P_sec] = calculate_Pcond_3Layer( ...
        params.I, params.f, ...
        Rdc_pri, Rdc_sec, ...
        params.t_cu_pri, params.t_cu_sec, ...
        params.rho_cu, params.u0);

    %% Total Loss
    P_total = P_core + P_pri + P_sec;

    %% Package Results
    result = struct( ...
        'Bmax',         Bmax,           ...
        'Ac',           Ac,             ...
        'r_centerpost', r_centerpost,   ...
        'w_core',       w_core,         ...
        'l_leg',        l_leg,          ...
        'h_yoke',       h_yoke,         ...
        'l_core',       l_core,         ...
        'h_core',       h_core,         ...
        'w_b',          params.w_b,     ...
        'l_winding',    l_winding,      ...
        'w_winding',    w_winding,      ...
        'V_total',      V_total,        ...
        'A_footprint',  A_footprint,    ...
        'P_core',       P_core,         ...
        'Rdc_pri',      Rdc_pri,        ...
        'Rdc_sec',      Rdc_sec,        ...
        'P_pri',        P_pri,          ...
        'P_sec',        P_sec,          ...
        'P_total',      P_total         ...
    );

end