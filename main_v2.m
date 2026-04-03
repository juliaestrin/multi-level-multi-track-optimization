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


%%  DESIGN SPECIFICATIONS
fprintf('\n--- Design Specifications ---\n');

% --- Electrical Specifications ---
Vin_nom     = 1500;         % [V]    Nominal input voltage
Vo_nom      = 48;           % [V]    Nominal output voltage per track
nt          = 2;            % [-]    Number of primary tracks (series connection)
Pmax        = 6.25e3;       % [W]    Maximum output power
Pmin        = 0.1 * Pmax;   % [W]    Minimum output power (10% load)

% --- Topology/Frequency Selection ---
% topology    = "Multitrack";
% fsw         = 1000e3;        % [Hz]   FCML switching frequency
% f0          = fsw;           % [Hz]   Transformer frequency
% SiCData     = 'SiC Data Multitrack.xlsx';
% GaNData     = [];

topology    = "Multilevel Multitrack";
fsw         = 500e3;        % [Hz]   FCML switching frequency
f0          = 2*fsw;        % [Hz]   Transformer frequency
SiCData     = 'SiC Data tf.xlsx';
GaNData     = 'GaN Data tf.xlsx';
                            
% --- LLC Resonant Tank Specifications ---
Mg_nom      = 1.0;          % [-]    Nominal LLC gain (unity at resonance)
percentReg  = 0.1;          % [-]    Line regulation tolerance (±10%)
f_per       = 0.25;         % [-]    Frequency range (±25%)
Ln          = 5;            % [-]    Inductance ratio Lm/Lr


% --- Transformer Design Specifications ---
material_name = 'F80';      % Core material selection
                            % Options: 'ML91S' (high freq), 'F80' (general),
                            %          'N87', 'N97', '3F4'

w_h_max       = 20e-3;      % [m]    Maximum window height constraint
w_scale       = 1;          % [-]    Winding width scale factor
                            %        1.0 = square winding (w = l)
                            %        0.5 = rectangular (w = 0.5*l)
h_core_max  = 100e-3;       % [m]
w_max = 200e-3;
l_max = 200e-3; 

s_ct = 0.5e-3;      % core to trace spacing
h_pcb = 1.6e-3;   % pcb height

centerpost_shape = 'round'; 
stackup       = '5layer' ;   % Winding layer configuration
%   Supported stackup configurations:
%     '3layer'             - 3-layer: P-S-P
%     '5layer'             - 5-layer: P-P-P-P-S
%     '5layer_interleaved' - 5-layer interleaved: P-P-S-P-P
%     '6layer'             - 6-layer non-interleaved: P-P-S-S-P-P
%     '6layer_interleaved' - 6-layer interleaved: P-S-P-P-S-P
%     '7layer_interleaved' - 7-layer interleaved: P-S-P-S-P-S-P
%     '8layer_interleaved' - 8-layer interleaved: P-S-P-S-P-S-P-S


% --- Transformer Fixed Parameters ---
t_cu_pri    = 2 * 35e-6;    % [m]    Primary copper thickness (2 oz)
t_cu_sec    = 2 * 35e-6;    % [m]    Secondary copper thickness (2 oz)


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


% fprintf('  Input voltage:    %.0f V\n', Vin_nom);
% fprintf('  Output voltage:   %.0f V (per track)\n', Vo_nom);
% fprintf('  Tracks (series):  %d\n', nt);
% fprintf('  Power (max):      %.2f kW\n', Pmax / 1e3);
% fprintf('  Power (min):      %.2f kW\n', Pmin / 1e3);
% fprintf('  fsw (FCML):       %.0f kHz\n', fsw / 1e3);
% fprintf('  f0 (transformer): %.1f MHz\n', f0 / 1e6);
% fprintf('  LLC Mg_nom:       %.1f\n', Mg_nom);
% fprintf('  Line reg:         ±%.0f%%\n', percentReg * 100);
% fprintf('  Freq range:       ±%.0f%%\n', f_per * 100);
% fprintf('  Ln (Lm/Lr):       %.1f\n', Ln);
% fprintf('  Core material:    %s\n', material_name);
% fprintf('  Max window height:%.1f mm\n', w_h_max * 1e3);
% fprintf('  Winding scale:    %.2f\n', w_scale);
% fprintf('  Stackup:          %s\n', stackup);
% fprintf('  Copper (pri):     %.0f µm (%.1f oz)\n', t_cu_pri * 1e6, t_cu_pri / 35e-6);
% fprintf('  Copper (sec):     %.0f µm (%.1f oz)\n', t_cu_sec * 1e6, t_cu_sec / 35e-6);

%%  LLC RESONANT TANK DESIGN

fprintf('\n --- LLC RESONANT TANK DESIGN ---\n');

LLC_design = designLLC_v2(topology, Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, ...
                       f_per, Pmax, Pmin, Ln);

% --- Extract LLC Design Results ---
Lu     = LLC_design.Lm;         % [H]  Magnetizing inductance
Llk    = LLC_design.Lr;         % [H]  Resonant (leakage) inductance
Ir_rms = LLC_design.Ir_rms;     % [A]  Primary RMS current (worst case)
Ir_pk  = Ir_rms * sqrt(2);      % [A]  Primary peak current
n      = LLC_design.N;          % [-]  Turns ratio (primary:secondary)
np     = n / nt;                % [-]  Track turns ratio 

% fprintf('  Turns ratio:      %d:1\n', n);
% fprintf('  Np:1/2:           %d turns\n', np);
% fprintf('  Lr (resonant):    %.4f µH\n', Llk * 1e6);
% fprintf('  Lm (magnetizing): %.4f µH\n', Lu * 1e6);
% fprintf('  Cr:               %.4f µF\n', LLC_design.Cr * 1e6);
% fprintf('  Ir_rms:           %.2f A\n', Ir_rms);
% fprintf('  Ir_pk:            %.2f A\n', Ir_pk);
% fprintf('  Qe_max:           %.4f\n', LLC_design.Qe_max);


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
Bsat  = material.Bsat; 

% --- Package Design Parameters ---
% Bundle all transformer design parameters into a struct
design_params = struct( ...
    'topology', topology,...
    'Vo',       Vo_nom,     ...  % Output voltage per track
    'nt',       nt,         ...  % Number of secondary tracks
    'Lu',       Lu,         ...  % Magnetizing inductance (from LLC design)
    'I',        Ir_pk,      ...  % Peak primary current
    'f',        f0,         ...  % Operating frequency (transformer core)
    'np',       np,         ...  % Primary turns ratio per track
    'k',        k,          ...  % Steinmetz coefficient
    'beta',     beta,       ...  % Steinmetz exponent
    'alpha',    alpha,      ...  % Steinmetz alpha (for modified Steinmetz)
    'uc',       uc,         ...  % Core relative permeability
    'Bsat',     Bsat,       ...
    'rho_cu',   rho_cu,     ...  % Copper resistivity
    'sigma_cu', sigma_cu,   ...  % Copper conductivity
    'u0',       u0,         ...  % Permeability of free space
    't_cu_pri', t_cu_pri,   ...  % Primary copper thickness
    't_cu_sec', t_cu_sec,   ...  % Secondary copper thickness
    'stackup',  stackup,     ...  % Winding layer configuration
    'centerpost_shape',  centerpost_shape,     ...  % Winding layer configuration
    's_ct', s_ct, ...
    'h_pcb', h_pcb ... 
    );


% --- Run Optimization ---

[x_opt, P_opt, results] = optimize_VIRT_v2(w_max, l_max, design_params)

T_tx = calculate_transformer_temp(results.P_total, results.Ac, results.h_w, R_plate, Area_plate, T_water, sig_grease, d_grease)


idx = 0; 
VIRT3Dfigure_ROUNDCenter_V2(results, material, T_tx, idx)
%VIRT3Dfigure_v2(results, material, T_tx, idx)

% % --- Display Optimal Design Results ---
% fprintf('\nOptimal Transformer Design:\n');
% fprintf('  Pv_max_opt:       %.2f kW/m³\n', opt.Pv_max_opt / 1e3);
% fprintf('  w_height_opt:     %.2f mm\n', opt.w_height_opt * 1e3);
% fprintf('  Bmax_opt:         %.4f T\n', opt.Bmax_opt);
% fprintf('  Core volume:      %.2f cm³\n', opt.opt_design.V_total * 1e6);
% fprintf('  Footprint:        %.2f cm²\n', opt.opt_design.A_footprint * 1e4);
% fprintf('  Air gap:          %.4f mm\n', opt.opt_design.lg * 1e3);
% fprintf('  P_core:           %.2f W\n', opt.P_core_min);
% fprintf('  P_copper:         %.2f W\n', opt.P_copper_min);
% fprintf('  P_total:          %.2f W\n', opt.P_total_min);
% fprintf('  T_transformer:    %.2f °C\n', T_tx);
% 
% % --- Generate 3D Visualization ---
% VIRT3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt, material.name, T_tx, 0 );
