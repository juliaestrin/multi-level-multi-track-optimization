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

% --- Heat Sink Parameters ---
R_plate = 0.08; % [C/W]
Area_plate = 152.4e-3 * 76.2e-3; % [m^2]
T_water = 45;


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
    fprintf('  Ir_rms:           %.4f A\n', Ir_rms);
    fprintf('  Ir_pk:            %.4f A\n', Ir_pk);
    fprintf('  Qe_max:           %.4f\n', LLC_design.Qe_max);

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
        Pmax, 0, 0);

    %% ---------- EXTRACT BEST POINT ----------
    [bestSummary, effTableAug] = extractBestPointSummary(effOut.table, topology);

    %% ---------- STORE RESULTS ----------
    results(ii).topology    = topology;
    results(ii).fsw         = fsw;
    results(ii).f0          = f0;
    results(ii).LLC         = LLC_design;
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