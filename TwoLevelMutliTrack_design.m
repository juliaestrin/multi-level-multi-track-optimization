% 2-Level Multitrack Converter Design, Optimization, and Pareto Front
%
% Authors: Julia Estrin, Qijia Li
%% 
%
% Description:
%   Runs the complete design workflow for the 2-Level Multitrack:
%   1. LLC resonant tank design
%   2. VIRT transformer optimization
%   3. Primary switch analysis
%   4. Secondary switch analysis
%   5. Overall efficiency / area pareto analysis

clear; close all; clc;

%% ========================= SETUP =========================
addpath('Switch Functions');
addpath('Transformer Design');
addpath('LLC Design');

%% ================= DESIGN SPECIFICATIONS =================

fprintf('\n--- Shared Design Specifications ---\n');

% --- Electrical Specifications ---
Vin_nom     = 1500;         % [V]    Nominal input voltage
Vo_nom      = 48;           % [V]    Nominal output voltage per track
nt          = 2;            % [-]    Number of primary tracks (series connection)
Pmax        = 6.25e3;       % [W]    Maximum output power
Pmin        = 0.1 * Pmax;   % [W]    Minimum output power (10% load)

% --- LLC Resonant Tank Specifications ---
Mg_nom      = 1.0;          % [-]    Nominal LLC gain
percentReg  = 0.05;          % [-]    Line regulation tolerance
f_per       = 0.25;         % [-]    Frequency range
Ln          = 5;            % [-]    Inductance ratio Lm/Lr

% --- Transformer Design Specifications ---
material_name = 'F80';      % Core material selection
                            % Options: 'ML91S' (high freq), 'F80' (general),
                            %          'N87', 'N97', '3F4'

w_max = 70e-3;             % [m] max transformer width [x-direction]
l_max = 70e-3;             % [m] max allowable transformer length [y-direction]
gap_loc = 'all';        % gap location: 'center' or 'all' legs

centerpost_shape = 'round';
stackup = '5Layer';
Tmax = 10000; 

% --- Transformer Fixed Parameters ---
t_cu_pri    = 2 * 35e-6;    % [m]    Primary copper thickness (2 oz)
t_cu_sec    = 2 * 35e-6;    % [m]    Secondary copper thickness (2 oz)
s_ct = 0.5e-3;              % [m]    core to trace spacing
h_pcb = 2.2e-3;             % [m]    pcb height

% --- Heat Sink Parameters ---
R_plate    = 0.08;                          % [C/W]
Area_plate = 152.4e-3 * 76.2e-3;           % [m^2]
T_water    = 45;                            % [C]

% --- Thermal Grease ---
sig_grease = 4;          % [W/(mK)]
d_grease   = 0.127e-3;   % [m]

% --- Material Constants ---
rho_cu      = 2.2e-8;        % [Ohm·m]
sigma_cu    = 1 / rho_cu;    % [S/m]
u0          = 4 * pi * 1e-7; % [H/m]

%% ================= TOPOLOGY CONFIGURATION =================

topology = "3-level Multitrack";
fsw      = 1000e3;              % [Hz]
f0       = fsw;                 % f0_factor = 1
SiCData  = 'SiC Data tf.xlsx';
GaNData  = 'GaN Data tf.xlsx';

fprintf('\n====================================================\n');
fprintf('TOPOLOGY: %s\n', topology);
fprintf('====================================================\n');
fprintf('  fsw (FCML):       %.0f kHz\n', fsw / 1e3);
fprintf('  f0 (transformer): %.1f MHz\n', f0 / 1e6);

%% ================= CORE MATERIAL =================

fprintf('\n--- Core Material ---\n');
material = get_core_material(material_name);

k     = material.k;
beta  = material.beta;
alpha = material.alpha;
uc    = material.uc;
Bsat  = material.Bsat;

fprintf('Core Material: %s\n', material.name);
fprintf('  k:                %.4e\n', k);
fprintf('  beta:             %.4f\n', beta);
if isfield(material, 'alpha')
    fprintf('  alpha:            %.4f\n', alpha);
end
fprintf('  uc:               %d\n', uc);

%% ---------- LLC RESONANT TANK DESIGN ----------
fprintf('\n--- LLC RESONANT TANK DESIGN ---\n');

LLC_design = designLLC_v3(topology, Vin_nom, Vo_nom, Mg_nom, nt, ...
    percentReg, fsw, f_per, Pmax, Pmin, Ln);

Lu     = LLC_design.Lm;
Llk    = LLC_design.Lr;
Ir_rms = LLC_design.Ir_rms;
Ir_pk  = Ir_rms * sqrt(2);
N      = LLC_design.N;
np     = N / nt;

fprintf('  Turns ratio:      %d:1\n', N);
fprintf('  Np:1/2:           %d turns\n', np);
fprintf('  Lr (resonant):    %.4f µH\n', Llk * 1e6);
fprintf('  Lm (magnetizing): %.4f µH\n', Lu * 1e6);
fprintf('  Cr:               %.4f µF\n', LLC_design.Cr * 1e6);
fprintf('  Ir_rms:           %.2f A\n', Ir_rms);
fprintf('  Ir_pk:            %.2f A\n', Ir_pk);
fprintf('  Qe_max:           %.4f\n', LLC_design.Qe_max);

%% ---------- PACKAGE DESIGN PARAMETERS ----------
design_params = struct( ...
    'topology',  topology,   ...
    'Vo',       Vo_nom,     ...
    'nt',       nt,         ...
    'Lu',       Lu,         ...
    'I',        Ir_pk,      ...
    'f',        f0,         ...
    'np',       np,         ...
    'k',        k,          ...
    'beta',     beta,       ...
    'alpha',    alpha,      ...
    'uc',       uc,         ...
    'Bsat',     Bsat,       ...
    'rho_cu',   rho_cu,     ...
    'sigma_cu', sigma_cu,   ...
    'u0',       u0,         ...
    't_cu_pri', t_cu_pri,   ...
    't_cu_sec', t_cu_sec,   ...
    'stackup',  stackup,    ...
    'centerpost_shape', centerpost_shape, ...
    's_ct',     s_ct,       ...
    'h_pcb',    h_pcb,       ...
    'T_max',    Tmax,        ... 
    'R_plate', R_plate, ...
    'Area_plate', Area_plate, ...
    'T_water', T_water, ...
    'sig_grease', sig_grease, ...
    'd_grease', d_grease, ...
    'gap_loc', gap_loc ...
    );

%% ---------- VIRT TRANSFORMER OPTIMIZATION ----------
fprintf('\n--- VIRT TRANSFORMER OPTIMIZATION ---\n');

[x_opt, P_opt, TX_design] = optimize_VIRT_2LMT(w_max, l_max, design_params);

T_tx = calculate_transformer_temp(TX_design.P_total, TX_design.Ac, TX_design.h_w, ...
    R_plate, Area_plate, T_water, sig_grease, d_grease);

Llk_TX = calc_Llk(design_params, TX_design); 
% Cps    = calc_Cps(TX_design);
% f_res  = 1 / (2 * pi * sqrt(Cps * Llk));

if strcmp(centerpost_shape, 'round')
    VIRT3Dfigure_round(TX_design, material, T_tx, 0);
else
    VIRT3Dfigure_square(TX_design, material, T_tx, 0);
end

%% ---------- PRIMARY SIDE SWITCH ANALYSIS ----------
fprintf('\n--- PRIMARY SIDE SWITCH ANALYSIS ---\n');

out1 = analyzePriSwitches_v3(topology, Pmax, fsw, Ir_rms, 1, 8, ...
    [], [], 10000, GaNData, SiCData);

try
    pareto   = out1.best.pareto;
    T_pareto = struct2table(pareto);

    wantedVars = {'device_name','jj','P_cond_W','P_off_W','P_gate_W','loss_W','area_mm2'};
    if all(ismember(wantedVars, T_pareto.Properties.VariableNames))
        T_pareto = T_pareto(:, wantedVars);
        T_pareto.Properties.VariableNames = { ...
            'DeviceName', 'ParallelCount', 'P_cond', 'P_off', ...
            'P_gate', 'TotalLoss', 'TotalArea_mm2'};
    end
    disp(T_pareto);
catch
    warning('Could not display primary pareto table.');
end

%% ---------- SECONDARY SIDE SWITCH ANALYSIS ----------
fprintf('\n--- SECONDARY SIDE SWITCH ANALYSIS ---\n');

out2 = analyzeSecSwitches(Pmax, f0, 1, 10, [4 6 8], [4 6 8]);

%% ---------- OVERALL EFFICIENCY / PARETO ----------
fprintf('\n--- OVERALL SYSTEM EFFICIENCY / PARETO ---\n');

effOut = calcEfficiency_v4(out1, out2, "pareto", "pareto", ...
    Pmax, TX_design.P_total, TX_design.A_footprint * 1e6);

[bestSummary, effTableAug] = extractBestPointSummary(effOut.table, topology);

%% ---------- RESULTS SUMMARY ----------

fprintf('\n====================================================\n');
fprintf('RESULTS SUMMARY\n');
fprintf('====================================================\n');

comparisonTable = table( ...
    string(topology), ...
    fsw/1e3, ...
    f0/1e3, ...
    LLC_design.N, ...
    LLC_design.Lr*1e6, ...
    LLC_design.Lm*1e6, ...
    LLC_design.Cr*1e6, ...
    LLC_design.Ir_rms, ...
    TX_design.Pv/1e3, ...
    TX_design.Bmax, ...
    TX_design.Vc*1e6, ...
    TX_design.A_footprint*1e6, ...
    TX_design.A_footprint*1e4, ...
    TX_design.A_footprint*(39.3701^2), ...
    TX_design.lg*1e3, ...
    TX_design.P_core, ...
    TX_design.P_pri, ...
    TX_design.P_sec, ...
    TX_design.P_total, ...
    T_tx, ...
    bestSummary.BestEfficiency_percent, ...
    bestSummary.BestSystemTotalLoss_W, ...
    bestSummary.BestSystemTotalArea_mm2, ...
    bestSummary.BestSystemTotalArea_in2, ...
    bestSummary.NumParetoPoints, ...
    'VariableNames', { ...
        'Topology', ...
        'fsw_kHz', ...
        'f0_kHz', ...
        'TurnsRatio_N', ...
        'Lr_uH', ...
        'Lm_uH', ...
        'Cr_uF', ...
        'Ir_rms_A', ...
        'PvmaxOpt_kWm3', ...
        'Bmax_T', ...
        'TransformerVolume_cm3', ...
        'TransformerFootprint_mm2', ...
        'TransformerFootprint_cm2', ...
        'TransformerFootprint_in2', ...
        'AirGap_mm', ...
        'TransformerCoreLoss_W', ...
        'TransformerPriCopperLoss_W', ...
        'TransformerSecCopperLoss_W', ...
        'TransformerTotalLoss_W', ...
        'TransformerTemp_C', ...
        'BestSystemEfficiency_percent', ...
        'BestSystemTotalLoss_W', ...
        'BestSystemTotalArea_mm2', ...
        'BestSystemTotalArea_in2', ...
        'NumParetoPoints' ...
    });

disp(comparisonTable);

fprintf('\n--- Transformer Summary ---\n');
fprintf('  Volume:      %.4f cm^3\n', TX_design.Vc * 1e6);
fprintf('  Temperature: %.2f C\n',    T_tx);
fprintf('  Core loss:   %.2f W\n',    TX_design.P_core);
fprintf('  Copper loss: %.2f W\n',    TX_design.P_pri + TX_design.P_sec);

fprintf('\n--- Augmented Pareto Table ---\n');
disp(effTableAug);

%% ================= PLOT: PARETO FRONT =================

figure('Name', 'Pareto Front: Area vs Loss');
hold on; grid on; box on;
plot(effTableAug.PlotArea_in2, effTableAug.PlotLoss_W, ...
    'LineStyle', '-', 'LineWidth', 1.8, 'Color', [0 0.4470 0.7410], ...
    'Marker', 'o', 'MarkerSize', 7, 'DisplayName', topology);
xlabel('Total Area [in$^2$]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Total Loss [W]',      'Interpreter', 'latex', 'FontSize', 15);
title('Pareto Front: Total Loss vs Total Area');
legend('Location', 'best');
set(gca, 'FontSize', 12);

figure('Name', 'Pareto Front: Area vs Efficiency');
hold on; grid on; box on;
plot(effTableAug.PlotArea_in2, effTableAug.PlotEfficiency_percent, ...
    'LineStyle', '-', 'LineWidth', 1.8, 'Color', [0 0.4470 0.7410], ...
    'Marker', 'o', 'MarkerSize', 7, 'DisplayName', topology);
xlabel('Total Area [in$^2$]', 'Interpreter', 'latex', 'FontSize', 15);
ylabel('Efficiency [\%]',     'Interpreter', 'latex', 'FontSize', 15);
title('Pareto Front: Efficiency vs Total Area');
legend('Location', 'best');
set(gca, 'FontSize', 12);

%% ================= OPTIONAL: SAVE TABLES =================

% writetable(comparisonTable, 'Topology_Comparison_Table.xlsx');
% writetable(effTableAug,     'ParetoTable_3LevelMultitrack.xlsx');

fprintf('\nDone.\n');

%% =========================================================
%% ================= LOCAL HELPER FUNCTIONS =================
%% =========================================================

function [bestSummary, tblOut] = extractBestPointSummary(tblIn, topologyName)

    tblOut   = tblIn;
    varNames = tblOut.Properties.VariableNames;

    effVar  = findExactOrDie(varNames, {'Efficiency'});
    areaVar = findExactOrDie(varNames, {'Area [mm2]', 'Area_mm2', 'Area_mm2_'});
    lossVar = findExactOrDie(varNames, {'Total Loss [W]', 'TotalLoss_W', 'Total_Loss_W', 'TotalLossW'});

    effData  = tblOut.(effVar);
    areaData = tblOut.(areaVar);
    lossData = tblOut.(lossVar);

    if max(effData) < 1.5
        effDataPercent = 100 * effData;
    else
        effDataPercent = effData;
    end

    mm2_to_in2 = 1 / (25.4^2);

    tblOut.PlotEfficiency_percent = effDataPercent;
    tblOut.PlotArea_mm2           = areaData;
    tblOut.PlotArea_in2           = areaData * mm2_to_in2;
    tblOut.PlotLoss_W             = lossData;

    [bestEff, bestIdx] = max(effDataPercent);

    bestSummary                         = struct();
    bestSummary.Topology                = string(topologyName);
    bestSummary.BestIdx                 = bestIdx;
    bestSummary.BestEfficiency_percent  = bestEff;
    bestSummary.BestSystemTotalArea_mm2 = areaData(bestIdx);
    bestSummary.BestSystemTotalArea_in2 = areaData(bestIdx) * mm2_to_in2;
    bestSummary.BestSystemTotalLoss_W   = lossData(bestIdx);
    bestSummary.NumParetoPoints         = height(tblOut);
    bestSummary.EfficiencyColumn        = string(effVar);
    bestSummary.AreaColumn              = string(areaVar);
    bestSummary.LossColumn              = string(lossVar);
end

function varName = findExactOrDie(varNames, candidates)
    varName = '';
    for i = 1:numel(candidates)
        hit = find(strcmp(varNames, candidates{i}), 1, 'first');
        if ~isempty(hit)
            varName = varNames{hit};
            return;
        end
    end
    error('Could not find any of these columns: %s', strjoin(candidates, ', '));
end