% The function to analyze the primary side switches

function out = analyzePriSwitches(Power, f_sw_typ, I_r, mode, max_para, selected_para, compare_list, max_total_loss, ganFile, sicFile)
% Supports:
%   - GaN-only  (sicFile empty)
%   - SiC-only  (ganFile empty)
%   - GaN+SiC   (both provided)
%
% Uses numeric column indexing, with separate column maps for GaN vs SiC.
% Plot style:
%   - GaN markers filled
%   - SiC markers empty
% Helper functions are minimally changed and remain compatible with your secondary-side code.

%% ===================== Defaults =====================
if nargin < 10, sicFile = 'SiC Data.xlsx'; end
if nargin < 9,  ganFile = 'GaN Data.xlsx'; end

if nargin < 8 || isempty(max_total_loss), max_total_loss = inf; end
if nargin < 7 || isempty(compare_list),  compare_list = [];  end

% Normalize empties
if isempty(ganFile), ganFile = []; end
if isempty(sicFile), sicFile = []; end

%% ===================== Decide jj_set =====================
if mode == 1
    jj_set = 1:max_para;
elseif mode == 2
    jj_set = selected_para(:).';
else
    error('mode must be 1 or 2');
end

%% ===================== Read tables (flexible: GaN-only / SiC-only / both) =====================
tbls = {};
techs = string.empty(1,0);

if ~isempty(ganFile)
    GaNData = readtable(ganFile);
    tbls{end+1} = GaNData;
    techs(end+1) = "GaN";
end

if ~isempty(sicFile)
    SiCData = readtable(sicFile);
    tbls{end+1} = SiCData;
    techs(end+1) = "SiC";
end

if isempty(tbls)
    error('You must provide at least one data file: ganFile or sicFile.');
end

% Total device count across all included technologies
n_sw = 0;
for t = 1:numel(tbls)
    n_sw = n_sw + height(tbls{t});
end

% Meta table (global index -> which table + row)
Meta = table('Size',[n_sw 4], ...
    'VariableTypes',["string","string","double","double"], ...
    'VariableNames',["Tech","Name","SourceID","RowInSource"]);

g = 0;
for sid = 1:numel(tbls)
    T = tbls{sid};
    tech = techs(sid);
    nT = height(T);

    for ii = 1:nT
        g = g + 1;
        Meta.Tech(g) = tech;
        Meta.SourceID(g) = sid;
        Meta.RowInSource(g) = ii;
        Meta.Name(g) = getDeviceName(T, ii, tech, g);
    end
end

%% ===================== Column mapping by index =====================
% NOTE: These indices are relative to each source table's column order.

% GaN indices
idx.GaN.Rth_jc = 4;
idx.GaN.Eloss  = 6;
idx.GaN.Coss   = 7;
idx.GaN.L_min  = 8;
idx.GaN.W_min  = 9;
idx.GaN.Rds    = 10;
idx.GaN.Tjmax  = 11;
idx.GaN.tf     = 13;
idx.GaN.Qg     = 14;
idx.GaN.Vg     = 15;

% SiC indices
idx.SiC.Rth_jc = 4;
idx.SiC.Eloss  = 10;
idx.SiC.Coss   = 11;
idx.SiC.L_min  = 12;
idx.SiC.W_min  = 13;
idx.SiC.Rds    = 14;
idx.SiC.Tjmax  = 15;
idx.SiC.tf     = 17;
idx.SiC.Qg     = 18;
idx.SiC.Vg     = 19;

% Safety: ensure mapping exists for chosen tech(s)
for k = 1:numel(techs)
    tech = techs(k);
    if ~isfield(idx, tech)
        error('No idx mapping defined for tech "%s".', tech);
    end
end

%% ===================== Known Parameters =====================
Vin = 1500; 
Vin_per = 0.1; 
Vin_max = Vin*(1+Vin_per);

Vout = 48;
R_via = 134.83;
spacing = 4*(0.00254);
D_via = 6*(0.00254);
T_water = 45;

R_plate = 0.08;
Area_plate = 152.4*76.2; % [mm2]

f_sw_max = f_sw_typ*(1.25);
f_sw_min = f_sw_typ*(0.75);

pri_sw_count = 8;

%% ===================== Current rating =====================
I_r_pk = I_r*sqrt(2);
I_rating = ceil(I_r_pk/5)*5;

%% ===================== Pre-allocate =====================
nJ   = numel(jj_set);

Area         = nan(n_sw,nJ);
P_cond       = nan(n_sw,nJ);
P_gate       = nan(n_sw,nJ);
P_off        = nan(n_sw,nJ);
P_total      = nan(n_sw,nJ);
P_total_plot = nan(n_sw,nJ);
T_c          = nan(n_sw,nJ);
T_j          = nan(n_sw,nJ);
T_j_plot     = nan(n_sw,nJ);

%% ===================== Main Loop =====================
for ii_global = 1:n_sw

    tech = Meta.Tech(ii_global);      % "GaN" or "SiC"
    sid  = Meta.SourceID(ii_global);  % which table in tbls
    ii   = Meta.RowInSource(ii_global);

    T = tbls{sid};
    map = idx.(tech);

    % ---- Read by numeric index ----
    L_min   = T{ii, map.L_min};
    W_min   = T{ii, map.W_min};

    Rth_jc   = T{ii, map.Rth_jc};
    R_ds_max = T{ii, map.Rds};
    %E_loss   = T{ii, map.Eloss};
    Coss     = T{ii, map.Coss};
    t_f      = T{ii, map.tf};
    Q_g      = T{ii, map.Qg};
    V_g      = T{ii, map.Vg};
    T_j_max  = T{ii, map.Tjmax};

    % ---- footprint-dependent via/board thermal ----
    N_L = floor(L_min*0.1/(D_via+2*spacing));
    N_W = floor(W_min*0.1/(D_via+2*spacing));
    N_vias_max = N_L*N_W;

    Area_vias = N_vias_max*pi*((D_via*10/2)^2);
    Area_fr4  = L_min*W_min - Area_vias;
    R_fr4     = 4350*1.6/Area_fr4;

    Rth_board_min = ((R_via/N_vias_max)*R_fr4/((R_via/N_vias_max)+R_fr4));

    for k = 1:nJ
        jj = jj_set(k);

        Area(ii_global,k) = pri_sw_count * jj * L_min * W_min;
        I_d = I_rating / jj;

        P_off(ii_global,k)  = (jj)*((I_d^2)*((t_f*(1e-9))^2)/(24*2*Coss*(1e-12)))/(1/(f_sw_max/2));

        duty_cond = 0.25;
        P_cond(ii_global,k) = duty_cond * jj * (I_d^2) * R_ds_max;
        P_gate(ii_global,k) = jj * V_g * Q_g * 1e-9 * (f_sw_max/2);

        P_total(ii_global,k) = P_off(ii_global,k) + P_cond(ii_global,k) + P_gate(ii_global,k);

        Rth_pw    = R_plate * Area_plate / (L_min * W_min);
        Rth_inter = 6*(0.6) / (L_min * W_min * 0.01);

        P_per_device = P_total(ii_global,k) / jj;

        T_c(ii_global,k) = P_per_device * (Rth_inter + Rth_board_min + Rth_pw) + T_water;
        T_j(ii_global,k) = P_per_device * (Rth_jc   + Rth_inter + Rth_board_min + Rth_pw) + T_water;

        thermal_check = T_j(ii_global,k) < (T_j_max - 30);
        loss_check    = pri_sw_count*P_total(ii_global,k) <= max_total_loss;

        if thermal_check && loss_check
            P_total_plot(ii_global,k) = pri_sw_count*P_total(ii_global,k);
            T_j_plot(ii_global,k)     = T_j(ii_global,k);
        else
            P_total_plot(ii_global,k) = NaN;
            T_j_plot(ii_global,k)     = NaN;
        end
    end
end

%% ===================== Colors =====================
baseColors = [
    0.1216    0.4667    0.7059
    1.0000    0.4980    0.0549
    0.1725    0.6275    0.1725
    0.8392    0.1529    0.1569
    0.5804    0.4039    0.7412
    0.5490    0.3373    0.2941
    0.8902    0.4667    0.7608
    0.4980    0.4980    0.4980
    0.7373    0.7412    0.1333
    0.0902    0.7451    0.8118
    0.7373    0.2235    0.2235
    0.2549    0.4118    0.8824
    0.4940    0.1840    0.5560
    0.9290    0.6940    0.1250
    0.3010    0.7450    0.9330
];

% Need enough colors for the largest RowInSource among all tech tables
maxIdx = max(Meta.RowInSource);

if size(baseColors,1) < maxIdx
    baseColors = lines(maxIdx);  % auto-expand palette if needed
else
    baseColors = baseColors(1:maxIdx,:);
end

% Per-device colors: indexed by RowInSource (resets per tech automatically)
colors = nan(n_sw,3);
for i = 1:n_sw
    colors(i,:) = baseColors(Meta.RowInSource(i),:);
end

%% ===================== jj -> marker mapping (ABSOLUTE) =====================
marker_by_jj = cell(1,100);
marker_by_jj{1}  = 'v';
marker_by_jj{2}  = '>';
marker_by_jj{3}  = 'd';
marker_by_jj{4}  = 'o';
marker_by_jj{5}  = '^';
marker_by_jj{6}  = 'p';
marker_by_jj{7}  = 's';
marker_by_jj{8}  = 'h';
marker_by_jj{9}  = '<';
marker_by_jj{10} = 'x';

%% ============================================================
% FIGURE 1: All devices  (GaN filled, SiC empty)
%% ============================================================
figure(1); clf;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
for i = 1:n_sw
    isGaN = (Meta.Tech(i) == "GaN");
    if isGaN
        ls = '-o';
        mfc = colors(i,:);
    else
        ls = ':s';
        mfc = 'none';
    end

    plot(jj_set, P_total_plot(i,:), ls, ...
        'Color', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'MarkerFaceColor', mfc, ...
        'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);

nexttile; hold on; grid on;
for i = 1:n_sw
    isGaN = (Meta.Tech(i) == "GaN");
    if isGaN
        ls = '-o';
        mfc = colors(i,:);
    else
        ls = ':s';
        mfc = 'none';
    end

    plot(jj_set, T_j_plot(i,:), ls, ...
        'Color', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'MarkerFaceColor', mfc, ...
        'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
ylabel('Junction Temperature [$^{\circ}C$]','Interpreter','latex','FontSize',15);
legend(Meta.Name,'Location','best','FontSize',10);

sgtitle(sprintf('Primary Devices (%s) ($P_{out}$=%.2f kW, $f_{sw}$=%.0f kHz)', ...
    strjoin(unique(Meta.Tech), " + "), Power/1e3, f_sw_typ/1e3), ...
    'Interpreter','latex','FontSize',15);

%% ============================================================
% FIGURE 2: BEST selection
%% ============================================================
rank.jj_rank = 6;
rank.N_keep  = 20;
rank.jj_list = [1 2 3 4 5 6 7 8];
rank.jj_list = resolveRankJjList(mode, rank.jj_list, max_para, selected_para);

C = buildSelectedTopN(Area, P_total_plot, jj_set, Meta, colors, rank);

plotLossAndAreaFigure_OnePower(C, rank, marker_by_jj, ...
    sprintf('Primary Devices TOP-%d (rank at %d-parallel)', C.N_total, rank.jj_rank), 2);

if isempty(compare_list)
    compare_list = jj_set; 
end

compare_list = resolveRankJjList(mode, compare_list, max_para, selected_para);

compare_list = compare_list(ismember(compare_list, jj_set));

if isempty(compare_list)
    compare_list = jj_set; 
end

Loss_Area_plot = Area .* P_total_plot;
best = pickBestDevices(Area, P_total_plot, Loss_Area_plot, jj_set, compare_list, Meta);

%% ===================== Outputs =====================
out = struct();
out.Power         = Power;
out.f_sw_typ      = f_sw_typ;
out.mode          = mode;
out.max_para      = max_para;
out.selected_para = selected_para;

out.Meta          = Meta;
out.jj_set        = jj_set;

out.Area          = Area;
out.P_total       = P_total;
out.T_j           = T_j;

out.P_total_plot  = P_total_plot;
out.T_j_plot      = T_j_plot;

out.rank          = rank;

out.compare_list_used = compare_list;
out.best = best;
out.max_total_loss = max_total_loss;

end