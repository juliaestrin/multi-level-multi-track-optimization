% Example calls:
% Primary Side: 
% 1) GaN only:
% out1 = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
%     [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
%     'GaN Data.xlsx', []);
%
% 2) SiC only:
% out1 = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
%     [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
%     [], 'SiC Data.xlsx');
%
% 3) Both:
out1 = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
    [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
    'GaN Data.xlsx', 'SiC Data.xlsx');
%
% Secondary Side
out2 = analyzeSecSwitches(6.25e3, 1000e3, 1, 10, [4 6 8], [4 6]);

%% Operating Point
Lu  = 0.63e-6 * 5;             % [H]  Inductance requirement for 6.25 kW
Llk = 0.63e-6; 
I   = 21.4 * sqrt(2);          % [A]  Peak primary current for 6.25 kW
fsw = 500e3;                   % [Hz] Switching frequency
f   = 2 * fsw;                 % [Hz] Transformer core frequency (2x switching due to FCML frequency doubling)
Np  = 4;                       % [-]  Number of primary turns (4 : 1/2 VIRT)
w_h_max = 100e-3;               % [m] max window height (width) 

%% Ferrite Material Parameters (Proterials ML91S)
k    = 1606769230.46174;        % Steinmetz coefficient k   (P = k * B^beta)
beta = 3.19564613808786;        % Steinmetz exponent beta
uc   = 900;                     % [-]  Relative permeability of core material

%% Copper Parameters
rho_cu   = 2.2e-8;              % [Ohm·m] Resistivity at 100°C
sigma_cu = 1 / rho_cu;          % [S/m]   Conductivity at 100°C
u0       = 4 * pi * 1e-7;       % [H/m]   Permeability of free space

stackup = '5layer'; 
%   Supported stackup configurations:
%     '3layer'             - 3-layer: P-S-P
%     '5layer'             - 5-layer: P-P-P-P-S
%     '5layer_interleaved' - 5-layer: P-P-S-P-P
%     '6layer'             - 6-layer non-interleaved: P-P-S-S-P-P
%     '6layer_interleaved' - 6-layer interleaved: P-S-P-P-S-P
%     '7layer_interleaved' - 7-layer interleaved: P-S-P-S-P-S-P
%     '8layer_interleaved' - 8-layer interleaved: P-S-P-S-P-S-P-S

%% Transformer Fixed Parameters
w_b      = 4e-3;                % [m]     Winding window breadth/depth
t_cu_pri = 2 * 35e-6;           % [m]     Primary copper thickness (2 oz)
t_cu_sec = 2 * 35e-6;           % [m]     Secondary copper thickness (2 oz)

%% Package Design Parameters
design_params = struct( ...
    'Lu',       Lu,       ...
    'I',        I,        ...
    'f',        f,        ...
    'Np',       Np,       ...
    'k',        k,        ...
    'beta',     beta,     ...
    'uc',       uc,       ...
    'rho_cu',   rho_cu,   ...
    'sigma_cu', sigma_cu, ...
    'u0',       u0,       ...
    'w_b',      w_b,      ...
    't_cu_pri', t_cu_pri, ...
    't_cu_sec', t_cu_sec,  ...
    'stackup', stackup ...
);

%% Optimization Sweep
Pv_max_list   = linspace(50e3, 500e3, 500);     % [W/m^3] core loss density sweep
w_height_list = linspace(5e-3, w_h_max, 50);     % [m]     window height sweep
opt = optimize_VIRT(Pv_max_list, w_height_list, design_params);

Cps = calculcate_Cps_3Layer(opt.opt_design.l_winding, opt.opt_design.w_winding, opt.opt_design.w_core); 
f_res = 1/(2*pi*sqrt(Cps * Llk));

% Efficiency Calculation
eff = calcEfficiency(out1, out2, "pareto", "minLoss", 6.25e3, opt.P_total_min, opt.opt_design.A_footprint*(10^6));

%% ===================== SWICTH ANALYSIS =====================

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

sec_sw_count = 4;

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

        Area(ii,k) = sec_sw_count * jj * L_min * W_min;

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
            P_total_plot(ii,k) = sec_sw_count * P_total(ii,k);
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
figure(3); clf;
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
    sprintf('Secondary side TOP-%d (rank at %d-parallel)', C.N_total, rank.jj_rank), 4);

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

%% ===================== EFFICIENCY CALCULATION =====================

function effOut = calcEfficiency(PriOut, SecOut, modePri, modeSec, Power, other_loss, other_area, marker_by_jj, figNum)
% Efficiency Calculation
    if nargin < 7 || isempty(other_loss), other_loss = 0; end
    if nargin < 8 || isempty(marker_by_jj), marker_by_jj = makeMarkerMap(); end
    if nargin < 9 || isempty(figNum), figNum = 10; end

    if ~isstring(modePri) || ~isscalar(modePri) || ~isstring(modeSec) || ~isscalar(modeSec)
        error('modePri and modeSec must be string scalars (e.g., "minLoss", "pareto").');
    end

    isParetoPri = (modePri == "pareto");
    isParetoSec = (modeSec == "pareto");

    if isParetoPri && isParetoSec
        error('Both modePri and modeSec cannot be "pareto" in this function.');
    end

    % ===== Case A: no pareto -> scalar efficiency only =====
    if ~isParetoPri && ~isParetoSec
        priOne = PriOut.best.(modePri);
        secOne = SecOut.best.(modeSec);

        loss_total = priOne.loss_W + secOne.loss_W + other_loss;
        area_total = priOne.area_mm2 + secOne.area_mm2 + other_area;

        eta = (Power - loss_total) / Power;

        effOut = struct();
        effOut.eta = eta;
        effOut.loss_total_W = loss_total;
        effOut.area_total_mm2 = area_total;
        effOut.pri = priOne;
        effOut.sec = secOne;
        return;
    end

    % ===== Case B: exactly one pareto -> compute list + plot =====
    if isParetoPri
        paretoSide = "pri";
        priCands = PriOut.best.pareto;
        secFixed = SecOut.best.(modeSec);
    else
        paretoSide = "sec";
        secCands = SecOut.best.pareto;
        priFixed = PriOut.best.(modePri);
    end

    if isParetoPri
        N = numel(priCands);
        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[],'jj',[],'pri',[],'sec',[]), N, 1);
        for k = 1:N
            loss_total = priCands(k).loss_W + secFixed.loss_W + other_loss;
            area_total = priCands(k).area_mm2 + secFixed.area_mm2 + other_area;
            eta = (Power - loss_total) / Power;

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj       = priCands(k).jj;   % marker encodes pareto-side jj
            pts(k).pri      = priCands(k);
            pts(k).sec      = secFixed;
        end
    else
        N = numel(secCands);
        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[],'jj',[],'pri',[],'sec',[]), N, 1);
        for k = 1:N
            loss_total = priFixed.loss_W + secCands(k).loss_W + other_loss;
            area_total = priFixed.area_mm2 + secCands(k).area_mm2;
            eta = (Power - loss_total) / Power;

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj       = secCands(k).jj;   % marker encodes pareto-side jj
            pts(k).pri      = priFixed;
            pts(k).sec      = secCands(k);
        end
    end

    effOut = struct();
    effOut.points = pts;
    effOut.paretoSide = paretoSide;

    % ===== Plot: Area vs Loss with jj markers + efficiency text =====
    figure(figNum); clf;
    hold on; grid on;

    msize = 90;
    N = numel(pts);
    colors = lines(N);

    devNames = strings(N,1); % Extract device names (pareto side device)
    
    for k = 1:N
        if paretoSide == "pri"
            devNames(k) = string(pts(k).pri.device_name);
        else
            devNames(k) = string(pts(k).sec.device_name);
        end
    end
    
    [uniqueDevs, ~, devIdx] = unique(devNames, 'stable');

    colors = lines(numel(uniqueDevs));

    % Plot
    for k = 1:N
        jj = pts(k).jj;

        mk = 'o';
        if jj <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj})
            mk = marker_by_jj{jj};
        end

        Color = colors(devIdx(k),:);

        scatter(pts(k).area_mm2, pts(k).loss_W, msize, ...
            'Marker', mk, 'MarkerFaceColor', Color, ...
            'MarkerEdgeColor', Color, 'LineWidth', 1.2);

        text(pts(k).area_mm2, pts(k).loss_W, ...
            sprintf('  %.2f%%', 100*pts(k).eta), ...
            'FontSize', 10, 'VerticalAlignment', 'middle');
    end

    xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);
    % title(sprintf('Area vs Loss (P=%.2f kW)', Power/1e3), ...
    %     'Interpreter','none');

    jjVals = unique([pts.jj]);
    h = gobjects(numel(jjVals),1);
    for i = 1:numel(jjVals)
        jj = jjVals(i);
        mk = 'o';
        if jj <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj})
            mk = marker_by_jj{jj};
        end
        h(i) = scatter(nan, nan, msize, 'Marker', mk, ...
            'MarkerFaceColor', 'k', ...
            'MarkerEdgeColor','k',...
            'LineWidth', 1.2);
    end
    legend(h, compose('%d-parallel', jjVals), ...
        'Location','best','Orientation','horizontal','FontSize',10);

end

%% ===================== LOCAL FUNCTIONS =====================

function marker_by_jj = makeMarkerMap()
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
end

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

    % ---- marker fill by tech if available ----
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