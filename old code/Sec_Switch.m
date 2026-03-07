out = analyzeSecSwitches(6.25e3, 1000e3, 1, 10, [4 6 8], [4 6]);

function out = analyzeSecSwitches(Power, f_sw_typ, mode, max_para, selected_para, compare_list, dataFile)
% Fig 1: ALL devices  -> (Loss vs jj) + (Tj vs jj)
% Fig 2: TOP-N only   -> (Loss vs jj) + (Area vs Loss)
%
% Also selects best designs restricted to jj in compare_list:
%   1) min Loss
%   2) min Area
%   3) min (Loss*Area)
%   4) Pareto front (Loss vs Area)
%
% Output 'best' ALWAYS includes both loss and area for each selection.

%% ===================== Defaults =====================
if nargin < 7 || isempty(dataFile)
    dataFile = 'Sec Data.xlsx';
end
if nargin < 6 || isempty(compare_list)
    compare_list = [];
end

%% ===================== Decide jj_set =====================
if mode == 1
    jj_set = 1:max_para;
elseif mode == 2
    jj_set = selected_para(:).';   % ensure row vector
else
    error('mode must be 1 or 2');
end

SecData = readtable(dataFile);

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
R_load = (Vout^2)/Power;

%% ===================== Current rating =====================
I_sec_pk = (pi/2)*(Vout/R_load);
I_rating = ceil(I_sec_pk/5)*5;
fprintf("the I_sec_pk is %d\n",I_sec_pk);
fprintf("the I_rating is %d\n",I_rating);

%% ===================== Pre-allocate =====================
n_sw = height(SecData);
nJ   = numel(jj_set);

Area         = nan(n_sw,nJ);
P_cond       = nan(n_sw,nJ);
P_gate       = nan(n_sw,nJ);
P_total      = nan(n_sw,nJ);
P_total_plot = nan(n_sw,nJ);
T_c          = nan(n_sw,nJ);
T_j          = nan(n_sw,nJ);
T_j_plot     = nan(n_sw,nJ);
T_per        = nan(n_sw,nJ);

%% ===================== Main Loop =====================
for ii = 1:n_sw
    L_min = SecData{ii,7};
    W_min = SecData{ii,8};

    N_L = floor(L_min*0.1/(D_via+2*spacing));
    N_W = floor(W_min*0.1/(D_via+2*spacing));
    N_vias_max = N_L*N_W;

    Area_vias = N_vias_max*pi*((D_via*10/2)^2);
    Area_fr4  = L_min*W_min - Area_vias;
    R_fr4     = 4350*1.6/Area_fr4;

    Rth_board_min = ((R_via/N_vias_max)*R_fr4/((R_via/N_vias_max)+R_fr4));

    Rth_jc   = SecData{ii,4};
    R_ds_max = SecData{ii,9};
    Q_g      = SecData{ii,13}; % [nC]
    V_g      = SecData{ii,14};
    T_j_max  = SecData{ii,10};

    for k = 1:nJ
        jj = jj_set(k);

        Area(ii,k) = jj * L_min * W_min;

        I_d = I_rating / jj;

        duty_cond = 0.5;
        P_cond(ii,k) = duty_cond * jj * (I_d^2) * R_ds_max;
        P_gate(ii,k) = jj * V_g * Q_g * 1e-9 * f_sw_max;

        P_total(ii,k) = P_cond(ii,k) + P_gate(ii,k);

        Rth_pw    = R_plate * Area_plate / (L_min * W_min);
        Rth_inter = 6*(0.6) / (L_min * W_min * 0.01);

        P_per_device = P_total(ii,k) / jj;

        T_c(ii,k) = P_per_device * (Rth_inter + Rth_board_min + Rth_pw) + T_water;
        T_j(ii,k) = P_per_device * (Rth_jc   + Rth_inter + Rth_board_min + Rth_pw) + T_water;

        if T_j(ii,k) < (T_j_max - 30)
            T_per(ii,k)        = (T_j_max - T_j(ii,k)) / T_j_max;
            P_total_plot(ii,k) = P_total(ii,k);
            T_j_plot(ii,k)     = T_j(ii,k);
        else
            P_total_plot(ii,k) = NaN;
            T_j_plot(ii,k)     = NaN;
        end
    end
end

%% ===================== Colors =====================
colors = [
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
if size(colors,1) < n_sw
    colors = lines(n_sw);
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
% FIGURE 1: All devices
%% ============================================================
figure(1); clf;
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
for i = 1:n_sw
    plot(jj_set, P_total_plot(i,:), '-o', ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);

nexttile; hold on; grid on;
for i = 1:n_sw
    plot(jj_set, T_j_plot(i,:), '-o', ...
        'Color', colors(i,:), ...
        'MarkerFaceColor', colors(i,:), ...
        'LineWidth', 1.5, 'MarkerSize', 6);
end
xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
ylabel('Junction Temperature [$^{\circ}C$]','Interpreter','latex','FontSize',15);
legend(SecData.Name,'Location','best','FontSize',10);

sgtitle(sprintf('Secondary side ($P_{out}$=%.2f kW, $f_{sw}$=%.0f kHz)', Power/1e3, f_sw_typ/1e3), ...
    'Interpreter','latex','FontSize',15);

%% ============================================================
% FIGURE 2: Top-N selection (Loss vs jj + Area vs Loss)
%% ============================================================

% ---- Top-N selection ----
rank.jj_rank = 8;         % ABSOLUTE jj used for ranking (must be in jj_set)
rank.N_keep  = 20;        % number of devices to keep
rank.jj_list = [4 6];   % markers used on RIGHT tile only

% Apply mode rule to jj_list (overlap behavior)
rank.jj_list = resolveRankJjList(mode, rank.jj_list, max_para, selected_para);

% Build selection (top-N by loss at jj_rank)
C = buildSelectedTopN(Area, P_total_plot, jj_set, SecData, colors, rank);

plotLossAndAreaFigure_OnePower(C, rank, marker_by_jj, ...
    sprintf('Secondary side TOP-%d (rank at %d-parallel)', C.N_total, rank.jj_rank), 2);

%% ============================================================
% Best Devices (min loss, min area, min loss*area, and Pareto)
%% ============================================================

if isempty(compare_list)
    compare_list = jj_set;
end

% Make compare_list compatible with mode
compare_list = resolveRankJjList(mode, compare_list, max_para, selected_para);

% Ensure compare_list values exist in jj_set
compare_list = compare_list(ismember(compare_list, jj_set));
if isempty(compare_list)
    compare_list = jj_set;
end

% Compute Loss*Area (product)
Loss_Area_plot = Area .* P_total_plot;

best = pickBestDevices(Area, P_total_plot, Loss_Area_plot, jj_set, compare_list, SecData);

%% ===================== Outputs =====================
out = struct();
out.Power         = Power;
out.f_sw_typ      = f_sw_typ;
out.mode          = mode;
out.max_para      = max_para;
out.selected_para = selected_para;

out.SecData       = SecData;
out.jj_set        = jj_set;

out.Area          = Area;
out.P_total       = P_total;
out.T_j           = T_j;

out.P_total_plot  = P_total_plot;
out.T_j_plot      = T_j_plot;
out.T_per         = T_per;

out.rank          = rank;

out.compare_list_used = compare_list;
out.best = best;

end

%% ============================================================
%                       LOCAL FUNCTIONS
%% ============================================================

function jj_list_out = resolveRankJjList(mode, jj_list_in, max_para, selected_para)
% mode=1: keep jj_list_in (clipped to 1:max_para)
% mode=2: overlap with selected_para if any; else fallback to selected_para

    jj_list_in = unique(jj_list_in(:).','stable');

    if mode == 1
        jj_list_out = jj_list_in(jj_list_in >= 1 & jj_list_in <= max_para);
        if isempty(jj_list_out)
            jj_list_out = 1:max_para;
        end

    elseif mode == 2
        sel = unique(selected_para(:).','stable');
        overlap = jj_list_in(ismember(jj_list_in, sel));

        if isempty(overlap)
            jj_list_out = sel;
        else
            jj_list_out = overlap;
        end

    else
        error('mode must be 1 or 2');
    end
end

function C = buildSelectedTopN(Area, P_total_plot, jj_set, SecData, colors, rank)
% Select top N_keep devices by P_total_plot(:, col_rank) ascending, ignoring NaNs.
% rank.jj_rank is ABSOLUTE jj value.

    col_rank = find(jj_set == rank.jj_rank, 1);
    if isempty(col_rank)
        error('rank.jj_rank=%d is not in jj_set=[%s].', rank.jj_rank, num2str(jj_set));
    end

    loss_ref = P_total_plot(:, col_rank);
    valid = ~isnan(loss_ref);

    [~, idx] = sort(loss_ref(valid), 'ascend');
    all_idx = find(valid);

    N_keep = min(rank.N_keep, numel(idx));
    sel = all_idx(idx(1:N_keep));

    C.Area       = Area(sel,:);
    C.P_loss     = P_total_plot(sel,:);
    C.PartNumber = string(SecData.Name(sel));
    C.N_total    = N_keep;
    C.jj_set     = jj_set;

    if size(colors,1) < numel(sel)
        C.colors = lines(numel(sel));
    else
        C.colors = colors(sel,:);
    end
end

function plotLossAndAreaFigure_OnePower(C, rank, marker_by_jj, figTitle, figNum)
% Figure: left = Loss vs jj (DEVICE legend)
%         right = Area vs Loss (jj marker legend ONLY)

    figure(figNum); clf;
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    jj_set = C.jj_set;

    jj_list = rank.jj_list;
    jj_list = jj_list(ismember(jj_list, jj_set));

    %% ---- Tile 1: Loss vs jj (TOP-N only) ----
    nexttile;
    hold on; grid on;

    for i = 1:C.N_total
        plot(jj_set, C.P_loss(i,:), ...
            'LineStyle', '-', ...
            'Marker', 'o', ...
            'Color', C.colors(i,:), ...
            'MarkerFaceColor', C.colors(i,:), ...
            'MarkerEdgeColor', C.colors(i,:), ...
            'LineWidth', 1.5, ...
            'MarkerSize', 6);
    end

    xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);
    legend(C.PartNumber, 'Location','best', 'FontSize',10);

    %% ---- Tile 2: Area vs Loss (jj legend only) ----
    nexttile;
    hold on; grid on;

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

            hh = scatter( ...
                C.Area(i,col), ...
                C.P_loss(i,col), ...
                msize, ...
                C.colors(i,:), ...
                mk, ...
                'filled', ...
                'MarkerEdgeColor', C.colors(i,:), ...
                'LineWidth', 1.2);

            if ~isgraphics(h_jj(k))
                h_jj(k) = hh;
            end
        end
    end

    xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);

    if ~isempty(jj_list)
        ok = isgraphics(h_jj);
        legend(h_jj(ok), ...
            compose('%d-parallel', jj_list(ok)), ...
            'Location','best', ...
            'Orientation','horizontal', ...
            'FontSize',10);
    end

    sgtitle(figTitle,'Interpreter','latex','FontSize',15);
end

function best = pickBestDevices(AreaFull, LossFull, LossAreaFull, jj_set, compare_list, SecData)
% Selects:
%  1) min Loss
%  2) min Area
%  3) min (Loss*Area)
%  4) Pareto front for (Area, Loss)
% Restricted to ABSOLUTE jj values in compare_list.
% Feasible points only (Loss not NaN).

    % Map compare_list ABSOLUTE jj -> column indices in full matrices
    cols = zeros(1, numel(compare_list));
    for k = 1:numel(compare_list)
        cols(k) = find(jj_set == compare_list(k), 1);
    end

    % Restricted-to-compare_list submatrices (still indexed by original device rows)
    Loss = LossFull(:, cols);          % n_sw x nCompare
    Area = AreaFull(:, cols);
    LossArea = LossAreaFull(:, cols);

    feas = ~isnan(Loss);               % feasible points

    % ---------- 1) Min Loss ----------
    Loss1 = Loss; Loss1(~feas) = NaN;
    [~, lin] = min(Loss1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Loss1), lin);
    jj_abs = compare_list(colLocal);
    best.minLoss = packOneLocal(dev, jj_abs, compare_list, Loss, Area, SecData);

    % ---------- 2) Min Area ----------
    Area1 = Area; Area1(~feas) = NaN;
    [~, lin] = min(Area1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Area1), lin);
    jj_abs = compare_list(colLocal);
    best.minArea = packOneLocal(dev, jj_abs, compare_list, Loss, Area, SecData);

    % ---------- 3) Min Loss*Area ----------
    LA1 = LossArea; LA1(~feas) = NaN;
    [~, lin] = min(LA1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(LA1), lin);
    jj_abs = compare_list(colLocal);
    best.minLossArea = packOneLocal(dev, jj_abs, compare_list, Loss, Area, SecData);

    % ---------- 4) Pareto front ----------
    [devIdx, colIdx] = find(feas);
    if isempty(devIdx)
        best.pareto = struct([]);   % empty struct array
        best.compare_list = compare_list;
        return;
    end

    area_vec = Area(sub2ind(size(Area), devIdx, colIdx));
    loss_vec = Loss(sub2ind(size(Loss), devIdx, colIdx));
    jj_vec   = compare_list(colIdx);

    pareto_mask = computeParetoFront(area_vec, loss_vec);

    % Indices of Pareto-optimal points in the (area_vec, loss_vec) list
    pDevIdx = devIdx(pareto_mask);
    pColIdx = colIdx(pareto_mask);

    % Build struct array: one struct per Pareto point
    paretoList = repmat(struct( ...
        'device_index', [], ...
        'device_name',  "", ...
        'jj',           [], ...
        'area_mm2',     [], ...
        'loss_W',       [], ...
        'loss_area_W_mm2_product', []), ...
        numel(pDevIdx), 1);

    for m = 1:numel(pDevIdx)
        di = pDevIdx(m);
        cj = pColIdx(m);

        jj_abs = compare_list(cj);
        a = Area(di, cj);
        l = Loss(di, cj);

        paretoList(m).device_index = di;
        paretoList(m).device_name  = string(SecData.Name(di));
        paretoList(m).jj           = jj_abs;
        paretoList(m).area_mm2     = a;
        paretoList(m).loss_W       = l;
        paretoList(m).loss_area_W_mm2_product = a * l;
    end

    % Optional: sort Pareto list by area (then loss) to make it nicer to read
    if ~isempty(paretoList)
        A = [paretoList.area_mm2].';
        L = [paretoList.loss_W].';
        [~, ord] = sortrows([A L], [1 2]);
        paretoList = paretoList(ord);
    end

    best.pareto = paretoList;
    best.compare_list = compare_list;
end

function s = packOneLocal(devIdx, jj_abs, compare_list, LossSub, AreaSub, SecData)
    colLocal = find(compare_list == jj_abs, 1);
    if isempty(colLocal)
        error('packOneLocal: jj_abs=%d not found in compare_list=[%s].', jj_abs, num2str(compare_list));
    end

    s = struct();
    s.device_index = devIdx;
    s.device_name  = string(SecData.Name(devIdx));
    s.jj           = jj_abs;

    s.loss_W  = LossSub(devIdx, colLocal);
    s.area_mm2 = AreaSub(devIdx, colLocal);
    s.loss_area_W_mm2_product = s.loss_W * s.area_mm2;
end

function isPareto = computeParetoFront(area_vec, loss_vec)
% Pareto-optimal points for minimizing (area, loss).
% i is dominated if exists j with area_j <= area_i AND loss_j <= loss_i
% and at least one strict.

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
