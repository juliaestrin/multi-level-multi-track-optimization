% Example calls:
% 1) GaN only:
% out = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
%     [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
%     'GaN Data.xlsx', []);
%
% 2) SiC only:
% out = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
%     [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
%     [], 'SiC Data.xlsx');
%
% 3) Both:
out = analyzePriSwitches(6.25e3, 800e3, 21.4142, 2, 10, ...
    [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
    'GaN Data.xlsx', 'SiC Data.xlsx');

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

        Area(ii_global,k) = jj * L_min * W_min;
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
        loss_check    = P_total(ii_global,k) <= max_total_loss;

        if thermal_check && loss_check
            P_total_plot(ii_global,k) = P_total(ii_global,k);
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

%% ===================== LOCAL FUNCTIONS =====================

function name = getDeviceName(T, ii, tech, fallbackIdx)
    if istable(T) && ~isempty(T) && ismember('Name', T.Properties.VariableNames)
        try
            nm = string(T.Name(ii));
            if strlength(nm) == 0 || nm == "missing"
                nm = "device#" + string(fallbackIdx);
            end
        catch
            nm = "device#" + string(fallbackIdx);
        end
    else
        nm = "device#" + string(fallbackIdx);
    end
    name = tech + ": " + nm;
end

function jj_list_out = resolveRankJjList(mode, jj_list_in, max_para, selected_para)
    jj_list_in = unique(jj_list_in(:).','stable');

    if mode == 1
        jj_list_out = jj_list_in(jj_list_in >= 1 & jj_list_in <= max_para);
        if isempty(jj_list_out), jj_list_out = 1:max_para; end
    elseif mode == 2
        sel = unique(selected_para(:).','stable');
        overlap = jj_list_in(ismember(jj_list_in, sel));
        if isempty(overlap), jj_list_out = sel; else, jj_list_out = overlap; end
    else
        error('mode must be 1 or 2');
    end
end

function C = buildSelectedTopN(Area, P_total_plot, jj_set, Meta, colors, rank)
    col_rank = find(jj_set == rank.jj_rank, 1);
    if isempty(col_rank)
        error('rank.jj_rank=%d is not in jj_set=[%s].', rank.jj_rank, num2str(jj_set));
    end

    loss_ref = P_total_plot(:, col_rank);
    valid = ~isnan(loss_ref);

    [~, idxs] = sort(loss_ref(valid), 'ascend');
    all_idx = find(valid);

    N_keep = min(rank.N_keep, numel(idxs));
    sel = all_idx(idxs(1:N_keep));

    C.Area       = Area(sel,:);
    C.P_loss     = P_total_plot(sel,:);
    C.PartNumber = string(Meta.Name(sel));
    C.N_total    = N_keep;
    C.jj_set     = jj_set;

    if size(colors,1) < numel(sel)
        C.colors = lines(numel(sel));
    else
        C.colors = colors(sel,:);
    end

    % ---- NEW (backward compatible): marker fill by tech if available ----
    if ismember('Tech', Meta.Properties.VariableNames)
        tech_sel = string(Meta.Tech(sel));
        C.MarkerFilled = (tech_sel == "GaN");   % GaN filled, SiC empty
    end
end

function plotLossAndAreaFigure_OnePower(C, rank, marker_by_jj, figTitle, figNum)
    figure(figNum); clf;
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    jj_set = C.jj_set;
    jj_list = rank.jj_list;
    jj_list = jj_list(ismember(jj_list, jj_set));

    nexttile; hold on; grid on;
    for i = 1:C.N_total

        fillFlag = true;
        if isfield(C,'MarkerFilled') && numel(C.MarkerFilled) >= i
            fillFlag = C.MarkerFilled(i);
        end

        if fillFlag
            ls = '-';
            mk = 'o';
            mfc = C.colors(i,:);
        else
            ls = ':';
            mk = 's';
            mfc = 'none';
        end

        plot(jj_set, C.P_loss(i,:), ...
            'LineStyle', ls, 'Marker', mk, ...
            'Color', C.colors(i,:), ...
            'MarkerFaceColor', mfc, ...
            'MarkerEdgeColor', C.colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    end
    xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);
    legend(C.PartNumber,'Location','best','FontSize',10);

    nexttile; hold on; grid on;
    msize = 90;
    h_jj = gobjects(numel(jj_list),1);

    for k = 1:numel(jj_list)
        jj_abs = jj_list(k);
        col = find(jj_set == jj_abs, 1);

        mk = 'o';
        if jj_abs <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_abs})
            mk = marker_by_jj{jj_abs};
        end

        for i = 1:C.N_total
            if isnan(C.P_loss(i,col)), continue; end

            fillFlag = true;
            if isfield(C,'MarkerFilled') && numel(C.MarkerFilled) >= i
                fillFlag = C.MarkerFilled(i);
            end

            if fillFlag
                hh = scatter(C.Area(i,col), C.P_loss(i,col), msize, ...
                    'Marker', mk, ...
                    'MarkerFaceColor', C.colors(i,:), ...
                    'MarkerEdgeColor', C.colors(i,:), ...
                    'LineWidth', 1.2);
            else
                hh = scatter(C.Area(i,col), C.P_loss(i,col), msize, ...
                    'Marker', mk, ...
                    'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', C.colors(i,:), ...
                    'LineWidth', 1.2);
            end

            if ~isgraphics(h_jj(k)), h_jj(k) = hh; end
        end
    end

    xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);

    if ~isempty(jj_list)
        ok = isgraphics(h_jj);
        legend(h_jj(ok), compose('%d-parallel', jj_list(ok)), ...
            'Location','best','Orientation','horizontal','FontSize',10);
    end

    sgtitle(figTitle,'Interpreter','latex','FontSize',15);
end

function best = pickBestDevices(AreaFull, LossFull, LossAreaFull, jj_set, compare_list, Meta)
    cols = zeros(1, numel(compare_list));
    for k = 1:numel(compare_list)
        cols(k) = find(jj_set == compare_list(k), 1);
    end

    Loss = LossFull(:, cols);
    Area = AreaFull(:, cols);
    LossArea = LossAreaFull(:, cols);
    feas = ~isnan(Loss);

    Loss1 = Loss; Loss1(~feas) = NaN;
    [~, lin] = min(Loss1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Loss1), lin);
    jj_abs = compare_list(colLocal);
    best.minLoss = packOneLocal(dev, jj_abs, compare_list, Loss, Area, Meta);

    Area1 = Area; Area1(~feas) = NaN;
    [~, lin] = min(Area1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Area1), lin);
    jj_abs = compare_list(colLocal);
    best.minArea = packOneLocal(dev, jj_abs, compare_list, Loss, Area, Meta);

    LA1 = LossArea; LA1(~feas) = NaN;
    [~, lin] = min(LA1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(LA1), lin);
    jj_abs = compare_list(colLocal);
    best.minLossArea = packOneLocal(dev, jj_abs, compare_list, Loss, Area, Meta);

    [devIdx, colIdx] = find(feas);
    if isempty(devIdx)
        best.pareto = struct([]);
        best.compare_list = compare_list;
        return;
    end

    area_vec = Area(sub2ind(size(Area), devIdx, colIdx));
    loss_vec = Loss(sub2ind(size(Loss), devIdx, colIdx));
    pareto_mask = computeParetoFront(area_vec, loss_vec);

    pDevIdx = devIdx(pareto_mask);
    pColIdx = colIdx(pareto_mask);

    paretoList = repmat(struct( ...
        'device_index', [], 'device_name', "", 'jj', [], ...
        'area_mm2', [], 'loss_W', [], 'loss_area_W_mm2_product', []), ...
        numel(pDevIdx), 1);

    for m = 1:numel(pDevIdx)
        di = pDevIdx(m);
        cj = pColIdx(m);

        jj_abs = compare_list(cj);
        a = Area(di, cj);
        l = Loss(di, cj);

        paretoList(m).device_index = di;
        paretoList(m).device_name  = string(Meta.Name(di));
        paretoList(m).jj           = jj_abs;
        paretoList(m).area_mm2     = a;
        paretoList(m).loss_W       = l;
        paretoList(m).loss_area_W_mm2_product = a * l;
    end

    if ~isempty(paretoList)
        A = [paretoList.area_mm2].';
        L = [paretoList.loss_W].';
        [~, ord] = sortrows([A L], [1 2]);
        paretoList = paretoList(ord);
    end

    best.pareto = paretoList;
    best.compare_list = compare_list;
end

function s = packOneLocal(devIdx, jj_abs, compare_list, LossSub, AreaSub, Meta)
    colLocal = find(compare_list == jj_abs, 1);
    if isempty(colLocal)
        error('packOneLocal: jj_abs=%d not found in compare_list=[%s].', jj_abs, num2str(compare_list));
    end

    s = struct();
    s.device_index = devIdx;
    s.device_name  = string(Meta.Name(devIdx));
    s.jj           = jj_abs;

    s.loss_W  = LossSub(devIdx, colLocal);
    s.area_mm2 = AreaSub(devIdx, colLocal);
    s.loss_area_W_mm2_product = s.loss_W * s.area_mm2;
end

function isPareto = computeParetoFront(area_vec, loss_vec)
    n = numel(area_vec);
    isPareto = true(n,1);
    for i = 1:n
        if ~isPareto(i), continue; end
        ai = area_vec(i);
        li = loss_vec(i);
        dominates_i = (area_vec <= ai) & (loss_vec <= li) & ((area_vec < ai) | (loss_vec < li));
        if any(dominates_i)
            isPareto(i) = false;
        end
    end
end