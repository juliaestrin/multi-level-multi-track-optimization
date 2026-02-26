% The function to analyze the secondary side switches

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