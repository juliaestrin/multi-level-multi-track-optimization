% Multilevel Multitrack Converter Design, Optimization, and Pareto Front
% Generation 
% Authors: Julia Estrin, Qijia Li
% Date:    03-02-2026
%
% Description:
%   Complete design workflow for the Multilevel Multitrack converter including 
%   front end Flying Capacitor Multilevel (FCML) and  VIRT transformer.
%   Performs LLC resonant tank design, VIRT transformer optimization, and 
%   semiconductor device selection for both primary and secondary sides.
%
% Design Flow:
%   1. Define specifications (voltage, power, frequency)
%   2. Design LLC resonant tank (Lr, Cr, Lm)
%   3. Optimize VIRT transformer (core geometry, losses)
%   4. Analyze primary side switches (GaN/SiC)
%   5. Analyze secondary side switches
%   6. Calculate overall system efficiency
%
% Key Features:
%   - FCML topology with frequency doubling (f_transformer = 2*f_switching)
%   - Multitrack (nt tracks in series)
%   - Configurable core material selection
%   - Configurable winding stackup (3-8 layers)
%   - 2D optimization sweep over Pv_max and window height
%   - Pareto front analysis for switch selection


clear; close all; clc;

fprintf('MULTILEVEL MULTITRACK CONVERTER DESIGN\n');

%%  SETUP - Add Required Function Paths
addpath('Switch Functions');     % Switch analysis functions
addpath('Transformer Design');   % VIRT transformer design functions
addpath('LLC Design');           % LLC resonant tank design functions


%  DESIGN SPECIFICATIONS
fprintf('\n--- Design Specifications ---\n');

% --- Electrical Specifications ---
Vin_nom     = 1500;         % [V]    Nominal input voltage
Vo_nom      = 48;           % [V]    Nominal output voltage per track
nt          = 2;            % [-]    Number of primary tracks (series connection)
Pmax        = 6.25e3;       % [W]    Maximum output power
Pmin        = 0.1 * Pmax;   % [W]    Minimum output power (10% load)


% --- Frequency Selection ---
topology    = "Multitrack";
fsw         = 1000e3;        % [Hz]   FCML switching frequency
f0          = fsw;           % [Hz]   Transformer frequency
SiCData     = 'SiC Data Multitrack.xlsx';
GaNData     = [];

% topology    = "Multilevel Multitrack";
% fsw         = 500e3;        % [Hz]   FCML switching frequency
% f0          = 2*fsw;        % [Hz]   Transformer frequency
% SiCData     = 'SiC Data tf.xlsx';
% GaNData     = 'GaN Data tf.xlsx';
                            

% --- LLC Resonant Tank Specifications ---
Mg_nom      = 1.0;          % [-]    Nominal LLC gain (unity at resonance)
percentReg  = 0.1;          % [-]    Line regulation tolerance (±10%)
f_per       = 0.25;         % [-]    Frequency range (±25%)
Ln          = 5;            % [-]    Inductance ratio Lm/Lr


% --- Transformer Design Specifications ---
material_name = 'F80';      % Core material selection
                            % Options: 'ML91S' (high freq), 'F80' (general),
                            %          'N87', 'N97', '3F4'

w_h_max       = 10e-3;      % [m]    Maximum window height constraint
w_scale       = 1;          % [-]    Winding width scale factor
                            %        1.0 = square winding (w = l)
                            %        0.5 = rectangular (w = 0.5*l)

stackup       = '5layer' ;   % Winding layer configuration
% Supported configurations:
%   '3layer'             - P-S-P
%   '5layer'             - P-P-P-P-S (non-interleaved)
%   '6layer'             - P-P-S-S-P-P (non-interleaved)
%   '6layer_interleaved' - P-S-P-P-S-P
%   '7layer_interleaved' - P-S-P-S-P-S-P
%   '8layer_interleaved' - P-S-P-S-P-S-P-S


% --- Transformer Fixed Parameters ---
w_b         = 4e-3;         % [m]    Winding window breadth (depth into page)
t_cu_pri    = 2 * 35e-6;    % [m]    Primary copper thickness (2 oz)
t_cu_sec    = 2 * 35e-6;    % [m]    Secondary copper thickness (2 oz)
h_core_max  = 24e-3;        % [m]

% --- Heat Sink Parameters ---
% PN: 180-10-6C
% https://wakefieldthermal.com/content/data_sheets/Standard%20Liquid%20Cold%20Plates.pdf
R_plate = 0.08; % [C/W] thermal resistance plate to inlet water
Area_plate = 152.4e-3*76.2e-3; % [m^2] area of cold plate 
T_water = 45; 

% Thermal Grease 
% BERGQUIST TGR 4000 (arbitrarily selected) 
sig_grease = 4;  % [W/(mK)] thermal conductivity of grease 
d_grease = 0.127e-3; % [m] reccomended thermal grease thickness 

% --- Material Constants ---
rho_cu      = 2.2e-8;       % [Ohm·m] Copper resistivity at 100°C
sigma_cu    = 1 / rho_cu;   % [S/m]   Copper conductivity
u0          = 4 * pi * 1e-7;% [H/m]   Permeability of free space


fprintf('  Input voltage:    %.0f V\n', Vin_nom);
fprintf('  Output voltage:   %.0f V (per track)\n', Vo_nom);
fprintf('  Tracks (series):  %d\n', nt);
fprintf('  Power (max):      %.2f kW\n', Pmax / 1e3);
fprintf('  Power (min):      %.2f kW\n', Pmin / 1e3);
fprintf('  fsw (FCML):       %.0f kHz\n', fsw / 1e3);
fprintf('  f0 (transformer): %.1f MHz\n', f0 / 1e6);
fprintf('  LLC Mg_nom:       %.1f\n', Mg_nom);
fprintf('  Line reg:         ±%.0f%%\n', percentReg * 100);
fprintf('  Freq range:       ±%.0f%%\n', f_per * 100);
fprintf('  Ln (Lm/Lr):       %.1f\n', Ln);
fprintf('  Core material:    %s\n', material_name);
fprintf('  Max window height:%.1f mm\n', w_h_max * 1e3);
fprintf('  Winding scale:    %.2f\n', w_scale);
fprintf('  Stackup:          %s\n', stackup);
fprintf('  Window breadth:   %.1f mm\n', w_b * 1e3);
fprintf('  Copper (pri):     %.0f µm (%.1f oz)\n', t_cu_pri * 1e6, t_cu_pri / 35e-6);
fprintf('  Copper (sec):     %.0f µm (%.1f oz)\n', t_cu_sec * 1e6, t_cu_sec / 35e-6);

%%  LLC RESONANT TANK DESIGN

fprintf('\n --- LLC RESONANT TANK DESIGN ---\n');


% Design LLC resonant tank based on gain requirements across
% the full line and load regulation range
% LLC_design = designLLC(Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, ...
%                        f_per, Pmax, Pmin, Ln);
LLC_design = designLLC_v2(topology, Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, ...
                       f_per, Pmax, Pmin, Ln);

% --- Extract LLC Design Results ---
Lu     = LLC_design.Lm;         % [H]  Magnetizing inductance
Llk    = LLC_design.Lr;         % [H]  Resonant (leakage) inductance
Ir_rms = LLC_design.Ir_rms;     % [A]  Primary RMS current (worst case)
Ir_pk  = Ir_rms * sqrt(2);      % [A]  Primary peak current
N      = LLC_design.N;          % [-]  Turns ratio (primary:secondary)
Np     = N / nt;                % [-]  Primary turns per track

fprintf('  Turns ratio:      %d:1\n', N);
fprintf('  Np:1/2:           %d turns\n', Np);
fprintf('  Lr (resonant):    %.4f µH\n', Llk * 1e6);
fprintf('  Lm (magnetizing): %.4f µH\n', Lu * 1e6);
fprintf('  Cr:               %.4f µF\n', LLC_design.Cr * 1e6);
fprintf('  Ir_rms:           %.2f A\n', Ir_rms);
fprintf('  Ir_pk:            %.2f A\n', Ir_pk);
fprintf('  Qe_max:           %.4f\n', LLC_design.Qe_max);


%%  VIRT TRANSFORMER OPTIMIZATION

fprintf('\n--- VIRT TRANSFORMER OPTIMIZATION ---\n');

% --- Load Core Material Properties ---
% Retrieves Steinmetz parameters (k, beta, alpha) and permeability (uc)
% for the selected ferrite material
material = get_core_material(material_name);

k     = material.k;             % Steinmetz coefficient k
beta  = material.beta;          % Steinmetz exponent beta
alpha = material.alpha;         % Steinmetz exponent alpha (if available)
uc    = material.uc;            % Relative permeability

fprintf('Core Material: %s\n', material.name);
fprintf('  k:                %.4e\n', k);
fprintf('  beta:             %.4f\n', beta);
if isfield(material, 'alpha')
    fprintf('  alpha:            %.4f\n', alpha);
end
fprintf('  uc:               %d\n', uc);

% --- Package Design Parameters ---
% Bundle all transformer design parameters into a struct
design_params = struct( ...
    'Vo',       Vo_nom,     ...  % Output voltage per track
    'nt',       nt,         ...  % Number of secondary tracks
    'Lu',       Lu,         ...  % Magnetizing inductance (from LLC design)
    'I',        Ir_pk,      ...  % Peak primary current
    'f',        f0,         ...  % Operating frequency (transformer core)
    'Np',       Np,         ...  % Primary turns per track
    'k',        k,          ...  % Steinmetz coefficient
    'beta',     beta,       ...  % Steinmetz exponent
    'alpha',    alpha,      ...  % Steinmetz alpha (for modified Steinmetz)
    'uc',       uc,         ...  % Core relative permeability
    'rho_cu',   rho_cu,     ...  % Copper resistivity
    'sigma_cu', sigma_cu,   ...  % Copper conductivity
    'u0',       u0,         ...  % Permeability of free space
    'w_b',      w_b,        ...  % Window breadth
    'w_scale',  w_scale,    ...  % Winding width scale factor
    't_cu_pri', t_cu_pri,   ...  % Primary copper thickness
    't_cu_sec', t_cu_sec,   ...  % Secondary copper thickness
    'stackup',  stackup     ...  % Winding layer configuration
);

% --- Define Optimization Sweep Ranges ---
% Sweep over core loss density (Pv_max) and window height (w_height)
% to find the design that minimizes total loss

Pv_max_list   = linspace(50e3, 1000e3, 1000-50+1);              % [W/m^3]
w_height_list = linspace(5e-3, w_h_max, (w_h_max*1e3-5)*2+1);  % [m]

fprintf('\nOptimization Sweep:\n');
fprintf('  Pv_max range:     %.0f - %.0f kW/m³ (%d points)\n', ...
    min(Pv_max_list)/1e3, max(Pv_max_list)/1e3, length(Pv_max_list));
fprintf('  w_height range:   %.1f - %.1f mm (%d points)\n', ...
    min(w_height_list)*1e3, max(w_height_list)*1e3, length(w_height_list));
fprintf('  Total designs:    %d\n', length(Pv_max_list) * length(w_height_list));

% --- Run Optimization ---
T_tx_max = 220;

opt = optimize_VIRT( ...
    Pv_max_list, w_height_list, ...
    h_core_max, T_tx_max, ...
    R_plate, Area_plate, T_water, sig_grease, d_grease, ...
    design_params);

% --- Approximate Transformer Temp Rise ---
T_tx = calculate_transformer_temp(opt.opt_design.P_total, opt.opt_design.Ac, opt.opt_design.h_core/2, R_plate, Area_plate, T_water, sig_grease, d_grease);

% --- Display Optimal Design Results ---
fprintf('\nOptimal Transformer Design:\n');
fprintf('  Pv_max_opt:       %.2f kW/m³\n', opt.Pv_max_opt / 1e3);
fprintf('  w_height_opt:     %.2f mm\n', opt.w_height_opt * 1e3);
fprintf('  Bmax_opt:         %.4f T\n', opt.Bmax_opt);
fprintf('  Core volume:      %.2f cm³\n', opt.opt_design.V_total * 1e6);
fprintf('  Footprint:        %.2f cm²\n', opt.opt_design.A_footprint * 1e4);
fprintf('  Air gap:          %.4f mm\n', opt.opt_design.lg * 1e3);
fprintf('  P_core:           %.2f W\n', opt.P_core_min);
fprintf('  P_copper:         %.2f W\n', opt.P_copper_min);
fprintf('  P_total:          %.2f W\n', opt.P_total_min);
fprintf('  T_transformer:    %.2f °C\n', T_tx);

% --- Generate 3D Visualization ---
core3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt, material.name, T_tx);

% --- Parasitic Capacitance and Resonance Check (Optional) ---
% Uncomment to calculate primary-to-secondary capacitance and verify
% that Llk-Cps resonance does not interfere with LLC operation
%
% Cps   = calculate_Cps_3Layer(opt.opt_design.l_winding, ...
%                               opt.opt_design.w_winding, ...
%                               opt.opt_design.w_core);
% f_res = 1 / (2 * pi * sqrt(Cps * Llk));
% %
% fprintf('  Cps:              %.4f nF\n', Cps * 1e9);
% fprintf('  f_res (Llk-Cps):  %.4f MHz\n', f_res * 1e-6);
% %
% if f_res < 1.5 * f0
%     warning('Resonance (%.2f MHz) too close to operating frequency (%.2f MHz)', ...
%         f_res / 1e6, f0 / 1e6);
% end


%%  PRIMARY SIDE SWITCH ANALYSIS (GaN and SiC)

% Evaluate different parallelization options (1-8 devices in parallel)
% for both GaN and SiC technologies
out1 = analyzePriSwitches_v2(topology, Pmax, fsw, Ir_rms, 1, 8, ...
    [], [1 2 3 4 5 6 7 8], 10000, ...
    GaNData, SiCData);

% Print out the loss breakdown for the MOSFETs on the pareto front
pareto = out1.best.pareto;

T_pareto = struct2table(pareto);

T_pareto = T_pareto(:, { ...
    'device_name', ...
    'jj', ...
    'P_cond_W', ...
    'P_off_W', ...
    'P_gate_W', ...
    'loss_W', ...
    'area_mm2'});

T_pareto.Properties.VariableNames = { ...
    'DeviceName', ...
    'ParallelCount', ...
    'P_cond', ...
    'P_off', ...
    'P_gate', ...
    'TotalLoss', ...
    'TotalArea_mm2'};

disp(T_pareto);

% out1 = analyzePriSwitches(Pmax, fsw, Ir_rms, 1, 8, ...
%     [], [1 2 3 4 5 6 7 8], 10000, ...
%     'GaN Data tf.xlsx', 'SiC Data tf.xlsx');

%%  SECONDARY SIDE SWITCH ANALYSIS

% Evaluate series-parallel combinations
% Series: 4, 6, 8 devices    Parallel: 4, 6, 8 devices
out2 = analyzeSecSwitches(Pmax, f0, 1, 10, [4 6 8], [4 6 8]);

% Calculate overall system efficiency and create pareto front including:
%   - Transformer losses/area (core + copper)
%   - Primary side switch losses/area
%   - Secondary side switch losses/area
%
% Device selection criteria:
%   modePri = "pareto"  → Evaluate all primary Pareto points
%   modeSec = "minLoss" → Use minimum loss secondary device

% Uncomment to run standard efficiency calculation:
% eff = calcEfficiency(out1, out2, "pareto", "minLoss", ...
%     Pmax, opt.P_total_min, opt.opt_design.A_footprint * 1e6);

% Run enhanced efficiency calculation with table output and Pareto analysis:
effOut = calcEfficiency_v3(out1, out2, "pareto", "pareto", ...
    Pmax, opt.P_total_min, opt.opt_design.A_footprint * 1e6);

% Sweep with switching frequency
% f_sw_list = [300e3, 400e3, 500e3];
% freqSweepOut = compareParetoFreq( ...
%     f_sw_list, ...
%     Pmax, Ir_rms, ...
%     Pv_max_list, w_height_list, design_params, ...
%     makeMarkerMap(), 20, false);

% Access results:
%   effOut.table          - Pareto front results table
%   effOut.points         - Individual design points
%   opt.opt_design        - Optimal transformer design struct
%   LLC_design            - LLC resonant tank design