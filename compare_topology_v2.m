% Multilevel Multitrack Converter Design, Optimization, and Pareto Front
% Generation - Multi-topology comparison version
%
% Authors: Julia Estrin, Qijia Li
% Modified for automated topology sweep + comparison
%
% Description:
%   Runs the complete design workflow for multiple topologies:
%   1. LLC resonant tank design
%   2. VIRT transformer optimization
%   3. Primary switch analysis
%   4. Secondary switch analysis
%   5. Overall efficiency / area pareto analysis
%   6. Cross-topology comparison table
%   7. Overlay pareto plots with distinct markers/colors

clear; close all; clc;

fprintf('MULTILEVEL MULTITRACK CONVERTER DESIGN\n');

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
percentReg  = 0.1;          % [-]    Line regulation tolerance
f_per       = 0.25;         % [-]    Frequency range
Ln          = 5;            % [-]    Inductance ratio Lm/Lr

% --- Transformer Design Specifications ---
material_name = 'F80';      % Core material selection
                            % Options: 'ML91S' (high freq), 'F80' (general),
                            %          'N87', 'N97', '3F4'

w_max = 50e-3;             % [m] max transformer width [x-direction]
w_scale       = 1;          % [-]    Winding width scale factor
                            %        1.0 = square winding (w = l)
                            %        0.5 = rectangular (w = 0.5*l)
l_max = 50e-3;             % [m] max allowable transformer length [y-direction]

centerpost_shape = 'round'; 
stackup       = '5layer' ;   % Winding layer configuration

% --- Transformer Fixed Parameters ---
t_cu_pri    = 2 * 35e-6;    % [m]    Primary copper thickness (2 oz)
t_cu_sec    = 2 * 35e-6;    % [m]    Secondary copper thickness (2 oz)
s_ct = 0.5e-3;              % [m]    core to trace spacing
h_pcb = 1.6e-3;             % [m]    pcb height

% --- Heat Sink Parameters ---
R_plate = 0.08; % [C/W]
Area_plate = 152.4e-3 * 76.2e-3; % [m^2]
T_water = 45;

% --- Thermal Grease ---
sig_grease = 4;          % [W/(mK)]
d_grease   = 0.127e-3;   % [m]

% --- Material Constants ---
rho_cu      = 2.2e-8;        % [Ohm·m]
sigma_cu    = 1 / rho_cu;    % [S/m]
u0          = 4 * pi * 1e-7; % [H/m]

T_tx_max = 150;

% fprintf('  Input voltage:      %.0f V\n', Vin_nom);
% fprintf('  Output voltage:     %.0f V (per track)\n', Vo_nom);
% fprintf('  Tracks (series):    %d\n', nt);
% fprintf('  Power (max):        %.2f kW\n', Pmax / 1e3);
% fprintf('  Power (min):        %.2f kW\n', Pmin / 1e3);
% fprintf('  LLC Mg_nom:         %.1f\n', Mg_nom);
% fprintf('  Line reg:           ±%.0f%%\n', percentReg * 100);
% fprintf('  Freq range:         ±%.0f%%\n', f_per * 100);
% fprintf('  Ln (Lm/Lr):         %.1f\n', Ln);
% fprintf('  Core material:      %s\n', material_name);
% fprintf('  Max window height:  %.1f mm\n', w_h_max * 1e3);
% fprintf('  Winding scale:      %.2f\n', w_scale);
% fprintf('  Stackup:            %s\n', stackup);
% fprintf('  Window breadth:     %.1f mm\n', w_b * 1e3);
% fprintf('  Copper (pri):       %.0f µm (%.1f oz)\n', t_cu_pri * 1e6, t_cu_pri / 35e-6);
% fprintf('  Copper (sec):       %.0f µm (%.1f oz)\n', t_cu_sec * 1e6, t_cu_sec / 35e-6);

%% ================= TOPOLOGY CONFIGURATION =================

topology_list = { ...
    % struct( ...
    %     "name", "Fullbridge LLC", ...
    %     "fsw", 1000e3, ...
    %     "f0_factor", 1, ...
    %     "SiCData", 'SiC Data Multitrack.xlsx', ...
    %     "GaNData", []), ...
    struct( ...
        "name", "Multitrack", ...
        "fsw", 1000e3, ...
        "f0_factor", 1, ...
        "SiCData", 'SiC Data Multitrack.xlsx', ...
        "GaNData", []), ...
    struct( ...
        "name", "3-level Multitrack", ...
        "fsw", 1000e3, ...
        "f0_factor", 1, ...
        "SiCData", 'SiC Data tf.xlsx', ...
        "GaNData", 'GaN Data tf.xlsx'), ...
    struct( ...
        "name", "Multilevel Multitrack", ...
        "fsw", 500e3, ...
        "f0_factor", 2, ...
        "SiCData", 'SiC Data tf.xlsx', ...
        "GaNData", 'GaN Data tf.xlsx') ...
    };

nTopo = numel(topology_list);

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

%% ================= OPTIMIZATION SWEEP =================

% Pv_max_list   = linspace(50e3, 1000e3, 1000 - 50 + 1);              % [W/m^3]
% w_height_list = linspace(5e-3, w_h_max, (w_h_max*1e3 - 5)*2 + 1);   % [m]
% 
% fprintf('\nOptimization Sweep:\n');
% fprintf('  Pv_max range:     %.0f - %.0f kW/m³ (%d points)\n', ...
%     min(Pv_max_list)/1e3, max(Pv_max_list)/1e3, length(Pv_max_list));
% fprintf('  w_height range:   %.1f - %.1f mm (%d points)\n', ...
%     min(w_height_list)*1e3, max(w_height_list)*1e3, length(w_height_list));
% fprintf('  Total designs:    %d\n', length(Pv_max_list) * length(w_height_list));
% 
%% ================= RESULTS PREALLOCATION =================

results = struct();

%% ================= MAIN TOPOLOGY LOOP =================

for ii = 1:nTopo

    cfg = topology_list{ii};

    topology = cfg.name;
    fsw      = cfg.fsw;
    f0       = cfg.f0_factor * fsw;
    SiCData  = cfg.SiCData;
    GaNData  = cfg.GaNData;

    fprintf('\n====================================================\n');
    fprintf('RUNNING TOPOLOGY: %s\n', topology);
    fprintf('====================================================\n');
    fprintf('  fsw (FCML):       %.0f kHz\n', fsw / 1e3);
    fprintf('  f0 (transformer): %.1f MHz\n', f0 / 1e6);

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
    % design_params = struct( ...
    %     'Vo',       Vo_nom,     ...
    %     'nt',       nt,         ...
    %     'Lu',       Lu,         ...
    %     'I',        Ir_pk,      ...
    %     'f',        f0,         ...
    %     'Np',       np,         ...
    %     'k',        k,          ...
    %     'beta',     beta,       ...
    %     'alpha',    alpha,      ...
    %     'uc',       uc,         ...
    %     'rho_cu',   rho_cu,     ...
    %     'sigma_cu', sigma_cu,   ...
    %     'u0',       u0,         ...
    %     'w_b',      w_b,        ...
    %     'w_scale',  w_scale,    ...
    %     't_cu_pri', t_cu_pri,   ...
    %     't_cu_sec', t_cu_sec,   ...
    %     'stackup',  stackup     ...
    % );
    design_params = struct( ...
        'topology',  topology,   ...
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
        'Bsat',     Bsat,       ...  % Core max B field
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

    %% ---------- VIRT TRANSFORMER OPTIMIZATION ----------
    fprintf('\n--- VIRT TRANSFORMER OPTIMIZATION ---\n');
    
    [x_opt, P_opt, TX_design] = optimize_VIRT(w_max, l_max, design_params);


    T_tx = calculate_transformer_temp(TX_design.P_total, TX_design.Ac, TX_design.h_w, R_plate, Area_plate, T_water, sig_grease, d_grease);

    % opt = optimize_VIRT( ...
    %     Pv_max_list, w_height_list, ...
    %     h_core_max, T_tx_max, ...
    %     R_plate, Area_plate, T_water, sig_grease, d_grease, ...
    %     design_params);
    % 
    % T_tx = calculate_transformer_temp( ...
    %     opt.opt_design.P_total, ...
    %     opt.opt_design.Ac, ...
    %     opt.opt_design.h_core/2, ...
    %     R_plate, Area_plate, T_water, sig_grease, d_grease);

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
    % try
    %     core3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt, ...
    %         sprintf('%s - %s', material.name, topology), T_tx);
    % catch ME
    %     warning('core3Dfigure failed for topology %s: %s', topology, ME.message);
    % end

    idx = ii; 
    if strcmp(centerpost_shape, 'round')
        VIRT3Dfigure_round(TX_design, material, T_tx, ii);
    else 
        VIRT3Dfigure_square(TX_design, material, T_tx, ii);
    end

    %% ---------- PRIMARY SIDE SWITCH ANALYSIS ----------
    fprintf('\n--- PRIMARY SIDE SWITCH ANALYSIS ---\n');

    out1 = analyzePriSwitches_v3(topology, Pmax, fsw, Ir_rms, 1, 8, ...
        [], [], 10000, GaNData, SiCData);

    try
        pareto = out1.best.pareto;
        T_pareto = struct2table(pareto);

        if topology == "Multilevel Multitrack"
            wantedVars = {'device_name','jj25','jj75','P_cond_W','P_off_W','P_gate_W','loss_W','area_mm2'};
            hasVars = ismember(wantedVars, T_pareto.Properties.VariableNames);

            if all(hasVars)
                T_pareto = T_pareto(:, wantedVars);
                T_pareto.Properties.VariableNames = { ...
                    'DeviceName', 'jj25', 'jj75', 'P_cond', 'P_off', ...
                    'P_gate', 'TotalLoss', 'TotalArea_mm2'};
            end
        else
            wantedVars = {'device_name','jj','P_cond_W','P_off_W','P_gate_W','loss_W','area_mm2'};
            hasVars = ismember(wantedVars, T_pareto.Properties.VariableNames);

            if all(hasVars)
                T_pareto = T_pareto(:, wantedVars);
                T_pareto.Properties.VariableNames = { ...
                    'DeviceName', 'ParallelCount', 'P_cond', 'P_off', ...
                    'P_gate', 'TotalLoss', 'TotalArea_mm2'};
            end
        end

        disp(T_pareto);
    catch
        warning('Could not display primary pareto table for %s.', topology);
    end

    %% ---------- SECONDARY SIDE SWITCH ANALYSIS ----------
    fprintf('\n--- SECONDARY SIDE SWITCH ANALYSIS ---\n');

    out2 = analyzeSecSwitches(Pmax, f0, 1, 10, [4 6 8], [4 6 8]);

    %% ---------- OVERALL EFFICIENCY / PARETO ----------
    fprintf('\n--- OVERALL SYSTEM EFFICIENCY / PARETO ---\n');

    effOut = calcEfficiency_v4(out1, out2, "pareto", "pareto", ...
        Pmax, TX_design.P_total, TX_design.A_footprint * 1e6);
    % effOut = calcEfficiency_v4(out1, out2, "pareto", "pareto", ...
    %     Pmax, 0, 0);

    %% ---------- EXTRACT BEST POINT ----------
    [bestSummary, effTableAug] = extractBestPointSummary(effOut.table, topology);

    %% ---------- STORE RESULTS ----------
    results(ii).topology    = topology;
    results(ii).fsw         = fsw;
    results(ii).f0          = f0;
    results(ii).LLC         = LLC_design;
    results(ii).TX_design   = TX_design;
    results(ii).T_tx        = T_tx;
    results(ii).out1        = out1;
    results(ii).out2        = out2;
    results(ii).effOut      = effOut;
    results(ii).effTableAug = effTableAug;
    results(ii).bestSummary = bestSummary;

end

%% ================= FULL COMPARISON TABLE =================

comparisonTable = table();

for ii = 1:nTopo
    r = results(ii);
    b = r.bestSummary;

        oneRow = table( ...
        string(r.topology), ...
        r.fsw/1e3, ...
        r.f0/1e3, ...
        r.LLC.N, ...
        r.LLC.Lr*1e6, ...
        r.LLC.Lm*1e6, ...
        r.LLC.Cr*1e6, ...
        r.LLC.Ir_rms, ...
        r.TX_design.Pv/1e3, ...
        r.TX_design.Bmax, ...
        r.TX_design.Vc*1e6, ...
        r.TX_design.A_footprint*1e6, ...
        r.TX_design.A_footprint*1e4, ...
        r.TX_design.A_footprint*(39.3701^2), ...
        r.TX_design.lg*1e3, ...
        r.TX_design.P_core, ...
        r.TX_design.P_pri, ...
        r.TX_design.P_sec, ...
        r.TX_design.P_total, ...
        r.T_tx, ...
        b.BestEfficiency_percent, ...
        b.BestSystemTotalLoss_W, ...
        b.BestSystemTotalArea_mm2, ...
        b.BestSystemTotalArea_in2, ...
        b.NumParetoPoints, ...
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

    comparisonTable = [comparisonTable; oneRow];
end

fprintf('\n====================================================\n');
fprintf('FULL COMPARISON TABLE\n');
fprintf('====================================================\n');
disp(comparisonTable);

%% ================= NEW TRANSFORMER COMPARISON TABLE =================
% Replaces the old "Best Point Comparison" figure

transformerComparisonTable = table();

for ii = 1:nTopo
    r = results(ii);

    oneRow = table( ...
        string(r.topology), ...
        r.TX_design.Vc * 1e6, ...
        r.T_tx, ...
        r.TX_design.P_core, ...
        r.TX_design.P_pri + r.TX_design.P_sec, ...
        'VariableNames', { ...
            'Topology', ...
            'TransformerVolume_cm3', ...
            'TransformerTemp_C', ...
            'P_core_W', ...
            'P_copper_W' ...
        });

    transformerComparisonTable = [transformerComparisonTable; oneRow];
end

fprintf('\n====================================================\n');
fprintf('TRANSFORMER COMPARISON TABLE\n');
fprintf('====================================================\n');
disp(transformerComparisonTable);

%% ================= DETAILED PARETO TABLES =================

for ii = 1:nTopo
    fprintf('\n----------------------------------------------------\n');
    fprintf('Augmented Pareto Table: %s\n', results(ii).topology);
    fprintf('----------------------------------------------------\n');
    disp(results(ii).effTableAug);
end

%% ================= PLOT 1: PARETO FRONT OVERLAY =================
% Area vs Loss, with area in in^2

figure('Name','Pareto Front Overlay: Area vs Loss');
hold on; grid on; box on;

styleMap = getTopologyPlotStyles();

for ii = 1:nTopo
    tbl = results(ii).effTableAug;
    sty = styleMap(ii);

    x = tbl.PlotArea_in2;
    y = tbl.PlotLoss_W;

    plot(x, y, ...
        'LineStyle', sty.LineStyle, ...
        'LineWidth', 1.8, ...
        'Color', sty.Color, ...
        'Marker', sty.Marker, ...
        'MarkerSize', 7, ...
        'DisplayName', results(ii).topology);
end

xlabel('Total Area [in$^2$]','Interpreter','latex','FontSize',15);
ylabel('Total Loss [W]','Interpreter','latex','FontSize',15);
title('Pareto Front Comparison: Total Loss vs Total Area');
legend('Location','best');
set(gca,'FontSize',12);

%% ================= PLOT 2: PARETO FRONT OVERLAY =================
% Area vs Efficiency, with area in in^2

figure('Name','Pareto Front Overlay: Area vs Efficiency');
hold on; grid on; box on;

for ii = 1:nTopo
    tbl = results(ii).effTableAug;
    sty = styleMap(ii);

    x = tbl.PlotArea_in2;
    y = tbl.PlotEfficiency_percent;

    plot(x, y, ...
        'LineStyle', sty.LineStyle, ...
        'LineWidth', 1.8, ...
        'Color', sty.Color, ...
        'Marker', sty.Marker, ...
        'MarkerSize', 7, ...
        'DisplayName', results(ii).topology);
end

xlabel('Total Area [in$^2$]','Interpreter','latex','FontSize',15);
ylabel('Efficiency [\%]','Interpreter','latex','FontSize',15);
title('Pareto Front Comparison: Efficiency vs Total Area');
legend('Location','best');
set(gca,'FontSize',12);

%% ================= OPTIONAL: SAVE TABLES =================

% writetable(comparisonTable, 'Topology_Comparison_Table.xlsx');
% writetable(transformerComparisonTable, 'Transformer_Comparison_Table.xlsx');

% for ii = 1:nTopo
%     fname = sprintf('ParetoTable_%s.xlsx', strrep(results(ii).topology, ' ', '_'));
%     writetable(results(ii).effTableAug, fname);
% end

fprintf('\nDone.\n');

%% =========================================================
%% ================= LOCAL HELPER FUNCTIONS =================
%% =========================================================

function styleMap = getTopologyPlotStyles()
    styleMap = struct([]);

    styleMap(1).Color     = [0 0.4470 0.7410];
    styleMap(1).Marker    = 'o';
    styleMap(1).LineStyle = '-';

    styleMap(2).Color     = [0.8500 0.3250 0.0980];
    styleMap(2).Marker    = 's';
    styleMap(2).LineStyle = '--';

    styleMap(3).Color     = [0.4660 0.6740 0.1880];
    styleMap(3).Marker    = '^';
    styleMap(3).LineStyle = '-.';

    styleMap(4).Color     = [0.4940 0.1840 0.5560];   % purple
    styleMap(4).Marker    = 'd';                      % diamond
    styleMap(4).LineStyle = ':';
end

function [bestSummary, tblOut] = extractBestPointSummary(tblIn, topologyName)

    tblOut = tblIn;
    varNames = tblOut.Properties.VariableNames;

    % Show actual variable names once if needed
    % disp(varNames')

    % ---- Find exact total columns first ----
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

    bestSummary = struct();
    bestSummary.Topology                = string(topologyName);
    bestSummary.BestIdx                 = bestIdx;
    bestSummary.BestEfficiency_percent  = bestEff;
    bestSummary.BestSystemTotalArea_mm2 = areaData(bestIdx);
    bestSummary.BestSystemTotalArea_in2 = areaData(bestIdx) * mm2_to_in2;
    bestSummary.BestSystemTotalLoss_W   = lossData(bestIdx);
    bestSummary.NumParetoPoints         = height(tblOut);

    bestSummary.EfficiencyColumn = string(effVar);
    bestSummary.AreaColumn       = string(areaVar);
    bestSummary.LossColumn       = string(lossVar);
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