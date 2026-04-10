% The function to analyze the primary side switches

function out = analyzePriSwitches_v3(topology, Power, f_sw_typ, I_r, mode, max_para, selected_para, compare_list, max_total_loss, ganFile, sicFile)
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
if nargin < 11, sicFile = 'SiC Data.xlsx'; end
if nargin < 10,  ganFile = 'GaN Data.xlsx'; end

if nargin < 9 || isempty(max_total_loss), max_total_loss = inf; end
if nargin < 8 || isempty(compare_list),  compare_list = [];  end

% Normalize empties
if isempty(ganFile), ganFile = []; end
if isempty(sicFile), sicFile = []; end

%% ===================== Decide candidate set =====================
if topology == "Multilevel Multitrack"
    % Each candidate is [jj25, jj75]
    if mode == 1
        % Full grid of all combinations from 1:max_para
        [J25, J75] = ndgrid(1:max_para, 1:max_para);
        jj_pair_set = [J25(:), J75(:)];
    elseif mode == 2
        % selected_para must be an N x 2 matrix: [jj25, jj75]
        if isempty(selected_para) || size(selected_para,2) ~= 2
            error(['For topology "Multilevel Multitrack", selected_para must be ', ...
                   'an N x 2 matrix with rows [jj25, jj75].']);
        end
        jj_pair_set = selected_para;
    else
        error('mode must be 1 or 2');
    end

    nJ = size(jj_pair_set,1);

else
    % Original scalar-jj behavior
    if mode == 1
        jj_set = 1:max_para;
    elseif mode == 2
        jj_set = selected_para(:).';
    else
        error('mode must be 1 or 2');
    end

    nJ = numel(jj_set);
end

%% ===================== Read tables (flexible: GaN-only / SiC-only / both) =====================
tbls = {};
techs = string.empty(1,0);

if ~isempty(ganFile)
    GaNData = readtable(ganFile,'VariableNamingRule','preserve');
    tbls{end+1} = GaNData;
    techs(end+1) = "GaN";
end

if ~isempty(sicFile)
    SiCData = readtable(sicFile,'VariableNamingRule','preserve');
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
idx.GaN.Cool   = 5;
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
idx.SiC.Cool   = 5;
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

if topology == "Multilevel Multitrack"
    pri_sw_count = 8;
    n25 = pri_sw_count/2;   % 4 positions at D = 0.25
    n75 = pri_sw_count/2;   % 4 positions at D = 0.75
elseif topology == "Multitrack"
    pri_sw_count = 4;
    duty_cond = 0.5;
elseif topology == "3-level Multitrack"
    pri_sw_count = 8;
    duty_cond = 0.5;
elseif topology == "Fullbridge LLC"
    pri_sw_count = 8;
    duty_cond = 0.5;
else
    error('Unsupported topology.');
end

%% ===================== Current rating =====================
I_r_pk = I_r*sqrt(2);
%I_rating = ceil(I_r_pk/5)*5;
I_rating = I_r_pk;

%% ===================== Pre-allocate =====================
Area         = nan(n_sw,nJ);
P_cond       = nan(n_sw,nJ);
P_gate       = nan(n_sw,nJ);
P_off        = nan(n_sw,nJ);
P_total      = nan(n_sw,nJ);
P_total_plot = nan(n_sw,nJ);
T_c          = nan(n_sw,nJ);
T_j          = nan(n_sw,nJ);
T_j_plot     = nan(n_sw,nJ);

% Optional but useful for multilevel case
T_j_25       = nan(n_sw,nJ);
T_j_75       = nan(n_sw,nJ);

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
    Coss     = T{ii, map.Coss};
    t_f      = T{ii, map.tf};
    Q_g      = T{ii, map.Qg};
    V_g      = T{ii, map.Vg};
    T_j_max  = T{ii, map.Tjmax};
    Cooling  = T{ii, map.Cool};

    % ---- footprint-dependent via/board thermal ----
    N_L = floor(L_min*0.1/(D_via+2*spacing));
    N_W = floor(W_min*0.1/(D_via+2*spacing));
    N_vias_max = N_L*N_W;

    Area_vias = N_vias_max*pi*((D_via*10/2)^2);
    Area_fr4  = L_min*W_min - Area_vias;
    R_fr4     = 4350*1.6/Area_fr4;

    if string(Cooling) == "Top"
        Rth_board_min = 0;
    else
        Rth_board_min = ((R_via/N_vias_max)*R_fr4/((R_via/N_vias_max)+R_fr4));
        %fprintf("The thermal resistance of PCB is %d\n", Rth_board_min);
    end 

    % Sweep
    for k = 1:nJ
        Rth_pw    = R_plate * Area_plate / (L_min * W_min);
        Rth_inter = 6*(0.6) / (L_min * W_min * 0.01);

        if topology == "Multilevel Multitrack"
            % Candidate pair: [jj25, jj75]
            jj25 = jj_pair_set(k,1);
            jj75 = jj_pair_set(k,2);

            I_d_25 = I_rating / jj25;
            I_d_75 = I_rating / jj75;

            % ----- Area -----
            Area(ii_global,k) = (n25*jj25 + n75*jj75) * L_min * W_min;

            % ----- Per-position loss for the D=0.25 group -----
            Poff_25  = jj25 * ((I_d_25^2) * ((t_f*1e-9)^2) / (24*2*Coss*1e-12)) / (1/f_sw_max);
            Pcond_25 = 0.25 * jj25 * (I_d_25^2) * R_ds_max;
            Pgate_25 = jj25 * V_g * Q_g * 1e-9 * f_sw_max;

            % ----- Per-position loss for the D=0.75 group -----
            Poff_75  = jj75 * ((I_d_75^2) * ((t_f*1e-9)^2) / (24*2*Coss*1e-12)) / (1/f_sw_max);
            Pcond_75 = 0.75 * jj75 * (I_d_75^2) * R_ds_max;
            Pgate_75 = jj75 * V_g * Q_g * 1e-9 * f_sw_max;

            % ----- Total over all switch positions -----
            P_off(ii_global,k)  = n25*Poff_25  + n75*Poff_75;
            P_cond(ii_global,k) = n25*Pcond_25 + n75*Pcond_75;
            P_gate(ii_global,k) = n25*Pgate_25 + n75*Pgate_75;

            P_total(ii_global,k) = P_off(ii_global,k) + P_cond(ii_global,k) + P_gate(ii_global,k);

            % ----- Thermal: evaluate worst group separately -----
            Pdev_25 = Poff_25 + Pcond_25 + Pgate_25;   % power per switch position of 0.25 group
            Pdev_75 = Poff_75 + Pcond_75 + Pgate_75;   % power per switch position of 0.75 group

            T_c_25 = Pdev_25 * (Rth_inter + Rth_board_min + Rth_pw) + T_water;
            T_c_75 = Pdev_75 * (Rth_inter + Rth_board_min + Rth_pw) + T_water;

            Tj_25 = Pdev_25 * (Rth_jc + Rth_inter + Rth_board_min + Rth_pw) + T_water;
            Tj_75 = Pdev_75 * (Rth_jc + Rth_inter + Rth_board_min + Rth_pw) + T_water;

            T_j_25(ii_global,k) = Tj_25;
            T_j_75(ii_global,k) = Tj_75;

            % Store worst-case for plotting
            T_c(ii_global,k) = max(T_c_25, T_c_75);
            T_j(ii_global,k) = max(Tj_25, Tj_75);

            thermal_check = (Tj_25 < (T_j_max - 30)) && (Tj_75 < (T_j_max - 30));
            loss_check    = P_total(ii_global,k) <= max_total_loss;

            if thermal_check && loss_check
                P_total_plot(ii_global,k) = P_total(ii_global,k);
                T_j_plot(ii_global,k)     = T_j(ii_global,k);
            else
                P_total_plot(ii_global,k) = NaN;
                T_j_plot(ii_global,k)     = NaN;
            end

        else
            % ===== Original scalar-jj logic for Multitrack =====
            jj = jj_set(k);
        
            Area(ii_global,k) = pri_sw_count * jj * L_min * W_min;
            I_d = I_rating / jj;

            P_off(ii_global,k)  = jj * ((I_d^2) * ((t_f*1e-9)^2) / (24*2*Coss*1e-12)) / (1/f_sw_max);
            P_cond(ii_global,k) = duty_cond * jj * (I_d^2) * R_ds_max;
            P_gate(ii_global,k) = jj * V_g * Q_g * 1e-9 * f_sw_max;

            P_total(ii_global,k) = P_off(ii_global,k) + P_cond(ii_global,k) + P_gate(ii_global,k);

            P_per_device = P_total(ii_global,k) / jj;

            T_c(ii_global,k) = P_per_device * (Rth_inter + Rth_board_min + Rth_pw) + T_water;
            T_j(ii_global,k) = P_per_device * (Rth_jc   + Rth_inter + Rth_board_min + Rth_pw) + T_water;

            thermal_check = T_j(ii_global,k) < (T_j_max - 30);
            loss_check    = pri_sw_count * P_total(ii_global,k) <= max_total_loss;

            if thermal_check && loss_check
                P_total_plot(ii_global,k) = pri_sw_count * P_total(ii_global,k);
                T_j_plot(ii_global,k)     = T_j(ii_global,k);
            else
                P_total_plot(ii_global,k) = NaN;
                T_j_plot(ii_global,k)     = NaN;
            end
        end
    end
end

%% STILL NEED TO FIX THE FOLLOWING PART!!!

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
if topology ~= "Multilevel Multitrack"
    %% ============================================================
    % FIGURE 1: All devices  (GaN filled, SiC empty)
    %% ============================================================
    figure(2); clf;
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
end

%% ============================================================
% FIGURE 2 / BEST selection
%% ============================================================
Loss_Area_plot = Area .* P_total_plot;

if topology == "Multilevel Multitrack"
    % ---------------------------------------------
    % skip buildSelectedTopN and plotting, because
    % they assume scalar jj rather than [jj25, jj75]
    % ---------------------------------------------
    rank = struct();
    rank.note = ['Skipped buildSelectedTopN / ranking plot for ', ...
                 'Multilevel Multitrack because this mode uses jj pairs [jj25, jj75].'];

    if isempty(compare_list)
        compare_list_used = jj_pair_set;
    else
        if size(compare_list,2) ~= 2
            error(['For topology "Multilevel Multitrack", compare_list must be ', ...
                   'an N x 2 matrix with rows [jj25, jj75].']);
        end

        tf_keep = ismember(compare_list, jj_pair_set, 'rows');
        compare_list_used = compare_list(tf_keep, :);

        if isempty(compare_list_used)
            compare_list_used = jj_pair_set;
        end
    end

    % In pair mode, P_cond / P_off / P_gate are already total losses
    best = pickBestDevices_v3( ...
        Area, P_total_plot, Loss_Area_plot, ...
        P_cond, P_off, P_gate, ...
        [], ...   % pri_sw_count not used in pair mode
        jj_pair_set, compare_list_used, Meta);

else
    % ---------------------------------------------
    % Original scalar-jj workflow
    % ---------------------------------------------
    rank.jj_rank = 6;
    rank.N_keep  = 20;
    rank.jj_list = [1 2 3 4 5 6 7 8];
    rank.jj_list = resolveRankJJList_v3(mode, rank.jj_list, max_para, selected_para);

    C = buildSelectedTopN(Area, P_total_plot, jj_set, Meta, colors, rank);

    plotLossAndAreaFigure_OnePower(C, rank, marker_by_jj, ...
        sprintf('Primary Devices TOP-%d (rank at %d-parallel)', C.N_total, rank.jj_rank), 3);

    if isempty(compare_list)
        compare_list_used = jj_set;
    else
        compare_list_used = resolveRankJJList_v2(mode, compare_list, max_para, selected_para);
        compare_list_used = compare_list_used(ismember(compare_list_used, jj_set));

        if isempty(compare_list_used)
            compare_list_used = jj_set;
        end
    end

    best = pickBestDevices_v3( ...
        Area, P_total_plot, Loss_Area_plot, ...
        P_cond, P_off, P_gate, ...
        pri_sw_count, ...
        jj_set, compare_list_used, Meta);
end
%best = pickBestDevices_v2(Area, P_total_plot, Loss_Area_plot, P_cond, P_off, P_gate, pri_sw_count, jj_set, compare_list, Meta);
%best = pickBestDevices(Area, P_total_plot, Loss_Area_plot, jj_set, compare_list, Meta);

%% ===================== Outputs =====================
out = struct();
out.Power         = Power;
out.f_sw_typ      = f_sw_typ;
out.mode          = mode;
out.max_para      = max_para;
out.selected_para = selected_para;

out.Meta          = Meta;
if topology == "Multilevel Multitrack"
    out.jj_pair_set = jj_pair_set;
    out.T_j_25 = T_j_25;
    out.T_j_75 = T_j_75;
else
    out.jj_set = jj_set;
end
out.Area          = Area;
out.P_total       = P_total;
out.T_j           = T_j;

out.P_total_plot  = P_total_plot;
out.T_j_plot      = T_j_plot;

out.rank          = rank;
out.topology      = topology;

out.compare_list_used = compare_list_used;
out.best = best;
out.max_total_loss = max_total_loss;

end