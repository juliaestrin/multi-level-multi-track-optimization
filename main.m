% =========================================================================
% MultiLevel Multitrack Converter Optimization 
% =========================================================================
% Authors: Julia Estrin, Qi Li
% Date:    03-02-2026
%
% DESIGN FLOW:
%   1. Define specifications → 2. Design LLC → 3. Optimize transformer →
%   4. Analyze switches → 5. Calculate efficiency
% =========================================================================

%% Add Required Paths
addpath('Switch Functions');
addpath('Transformer Design');
addpath('LLC Design');

%%  DESIGN SPECIFICATIONS

% --- Electrical Specifications ---
Vin_nom     = 1500;         % [V]    Nominal input voltage
Vo_nom      = 48;           % [V]    Nominal output voltage  
nt          = 2;            % [-]    Number of secondary tracks
Pmax        = 6.25e3;       % [W]    Maximum output power
Pmin        = 0.1 * Pmax;   % [W]    Minimum output power (10% load)

% --- Frequency Selection ---
fsw         = 400e3;        % [Hz]   FCML switching frequency
f0          = 2 * fsw;      % [Hz]   Transformer frequency (FCML doubles fsw)

% --- LLC Resonant Tank Specs ---
Mg_nom      = 1.0;          % [-]    Nominal LLC gain
percentReg  = 0.1;          % [-]    Line regulation ±10%
f_per       = 0.25;         % [-]    Frequency range ±25%
Ln          = 5;            % [-]    Inductance ratio Lm/Lr

% --- Transformer Design Specs ---
material_name = 'F80';      % Core: 'ML91S', 'F80', '3F46'
w_h_max       = 20e-3;     % [m]    Max window height
w_scale       = 1/2;        % [-]    Winding width scale (0.5 = half-width)
stackup       = '5layer';   % '3layer','5layer','6layer_interleaved', etc.

w_b         = 4e-3;         % [m]    Window breadth
t_cu_pri    = 2 * 35e-6;    % [m]    Primary copper (2 oz)
t_cu_sec    = 2 * 35e-6;    % [m]    Secondary copper (2 oz)

% --- Constants ---
rho_cu      = 2.2e-8;       % [Ohm·m] Copper resistivity at 100°C
sigma_cu    = 1 / rho_cu;   % [S/m]   Copper conductivity
u0          = 4 * pi * 1e-7;% [H/m]   Permeability of free space

%%  LLC RESONANT TANK DESIGN

LLC_design = designLLC(Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, ...
                       f_per, Pmax, Pmin, Ln);

% Extract results
Lu     = LLC_design.Lm;         % [H]  Magnetizing inductance
Llk    = LLC_design.Lr;         % [H]  Resonant inductance  
Ir_rms = LLC_design.Ir_rms;     % [A]  Primary RMS current
Ir_pk  = Ir_rms * sqrt(2);      % [A]  Primary peak current
N      = LLC_design.N;          % [-]  Turns Ratio 
Np     = N / nt;                % [-]  Total primary turns


%% VIRT TRANSFORMER OPTIMIZATION

% Load core material properties
material = get_core_material(material_name);
k    = material.k;              % Steinmetz coefficient
beta = material.beta;           % Steinmetz exponent
alpha = material.alpha;           % Steinmetz exponent
uc   = material.uc;             % Relative permeability

% Package parameters
design_params = struct( ...
    'Vo',       Vo_nom,     ...
    'nt',       nt,         ...
    'Lu',       Lu,         ...
    'I',        Ir_pk,      ...
    'f',        f0,         ...
    'Np',       Np,         ...
    'k',        k,          ...
    'beta',     beta,       ...
    'alpha',    alpha,       ...
    'uc',       uc,         ...
    'rho_cu',   rho_cu,     ...
    'sigma_cu', sigma_cu,   ...
    'u0',       u0,         ...
    'w_b',      w_b,        ...
    'w_scale',  w_scale,    ...
    't_cu_pri', t_cu_pri,   ...
    't_cu_sec', t_cu_sec,   ...
    'stackup',  stackup     ...
);

% Define sweep ranges
Pv_max_list   = linspace(50e3, 500e3, 500);     % [W/m^3]
w_height_list = linspace(5e-3, w_h_max, 50);    % [m]

% Run optimization
opt = optimize_VIRT(Pv_max_list, w_height_list, design_params);

% Visualize result
core3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt);

% Optional: Check parasitic resonance
% Cps   = calculate_Cps_3Layer(opt.opt_design.l_winding, ...
%                               opt.opt_design.w_winding, opt.opt_design.w_core);
% f_res = 1 / (2 * pi * sqrt(Cps * Llk));

%% PRIMARY SIDE SWITCH ANALYSIS (GaN and SiC)
out1 = analyzePriSwitches(Pmax, fsw, Ir_rms, 2, 10, ...
    [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
    'GaN Data.xlsx', 'SiC Data.xlsx');

%% SECONDARY SIDE SWITCH ANALYSIS
out2 = analyzeSecSwitches(Pmax, f0, 1, 10, [4 6 8], [4 6 8]);

%% SYSTEM EFFICIENCY CALCULATION
% eff = calcEfficiency(out1, out2, "pareto", "minLoss", ...
%     Pmax, opt.P_total_min, opt.opt_design.A_footprint * 1e6);

effOut = calcEfficiency_v2(out1, out2, "pareto", "minLoss", Pmax, opt.P_total_min, opt.opt_design.A_footprint * 1e6, [], [], 1);