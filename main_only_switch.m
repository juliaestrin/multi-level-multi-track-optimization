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

topology    = "2-level Multitrack";
fsw         = 1000e3;        % [Hz]   FCML switching frequency
f0          = fsw;        % [Hz]   Transformer frequency
SiCData     = 'SiC Data tf.xlsx';
GaNData     = 'GaN Data tf.xlsx';

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

% --- Heat Sink Parameters ---
% PN: 180-10-6C - https://wakefieldthermal.com/content/data_sheets/Standard%20Liquid%20Cold%20Plates.pdf
R_plate = 0.08; % [C/W] thermal resistance plate to inlet water
Area_plate = 152.4e-3*76.2e-3; % [m^2] area of cold plate 
T_water = 45; 

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

%%  LLC RESONANT TANK DESIGN

fprintf('\n --- LLC RESONANT TANK DESIGN ---\n');


% Design LLC resonant tank based on gain requirements across
% the full line and load regulation range
LLC_design = designLLC_v3(topology, Vin_nom, Vo_nom, Mg_nom, nt, ...
        percentReg, fsw, f_per, Pmax, Pmin, Ln);

% --- Extract LLC Design Results ---
Lu     = LLC_design.Lm;         % [H]  Magnetizing inductance
Llk    = LLC_design.Lr;         % [H]  Resonant (leakage) inductance
Ir_rms = LLC_design.Ir_rms;     % [A]  Primary RMS current (worst case)
Ir_pk  = Ir_rms * sqrt(2);      % [A]  Primary peak current
n      = LLC_design.N;          % [-]  Turns ratio (primary:secondary)
np     = n / nt;                % [-]  Track turns ratio 

fprintf('  Turns ratio:      %d:1\n', n);
fprintf('  Np:1/2:           %d turns\n', np);
fprintf('  Lr (resonant):    %.4f µH\n', Llk * 1e6);
fprintf('  Lm (magnetizing): %.4f µH\n', Lu * 1e6);
fprintf('  Cr:               %.4f µF\n', LLC_design.Cr * 1e6);
fprintf('  Ir_rms:           %.6f A\n', Ir_rms);
fprintf('  Ir_pk:            %.2f A\n', Ir_pk);
fprintf('  Qe_max:           %.4f\n', LLC_design.Qe_max);

%%  PRIMARY SIDE SWITCH ANALYSIS (GaN and SiC)

% Evaluate different parallelization options (1-8 devices in parallel)
% for both GaN and SiC technologies
out1 = analyzePriSwitches_v4(topology, Pmax, fsw, Ir_rms, 1, 8, ...
    [], [], 10000, ...
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
effOut = calcEfficiency_v4(out1, out2, "pareto", "pareto", ...
    Pmax, 0, 0);

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
%%

% design_summary_figure(topology, Vin_nom, Vo_nom, Pmax, fsw, f0, nt, ...
%        LLC_design, opt, material, stackup, t_cu_pri, t_cu_sec, T_tx);