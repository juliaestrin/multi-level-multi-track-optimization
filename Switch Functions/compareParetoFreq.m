function freqSweepOut = compareParetoFreq( ...
    f_sw_list, ...
    Pmax, Ir_rms, ...
    Pv_max_list, w_height_list, design_params, ...
    marker_by_jj, figNum, showEachTable)

% compareEfficiencyParetoAcrossFreq
%
% Sweep switching frequency, generate effOut at each frequency, and compare
% the Pareto-front points across frequencies.
%
% Plot rule for the comparison plot:
%   marker = primary parallel count
%   color  = switching frequency
%
% Inputs
%   f_sw_list      : vector of switching frequencies
%   Pmax           : output power
%   Ir_rms         : RMS ripple/current input for primary-side analysis
%   Pv_max_list    : optimizer sweep list
%   w_height_list  : optimizer sweep list
%   design_params  : struct for optimize_VIRT
%   marker_by_jj   : marker map cell array, optional
%   figNum         : comparison figure number, optional
%   showEachTable  : whether calcEfficiency_v3 shows table/figures, optional
%
% Output
%   freqSweepOut : struct containing all per-frequency results and the
%                  combined Pareto-front points used in the comparison plot

    % ---------------- Defaults ----------------
    if nargin < 7 || isempty(marker_by_jj)
        marker_by_jj = makeMarkerMap();
    end
    if nargin < 8 || isempty(figNum)
        figNum = 100;
    end
    if nargin < 9 || isempty(showEachTable)
        showEachTable = false;
    end

    nF = numel(f_sw_list);

    % Store each frequency result
    results = repmat(struct( ...
        'fsw', [], ...
        'out1', [], ...
        'out2', [], ...
        'opt', [], ...
        'effOut', []), nF, 1);

    % Combined points for plotting
    allPts = struct( ...
        'fsw', {}, ...
        'area_mm2', {}, ...
        'loss_W', {}, ...
        'eta', {}, ...
        'jj_pri', {}, ...
        'jj_sec', {}, ...
        'pri', {}, ...
        'sec', {});

    % ---------------- Sweep frequency ----------------
    for i = 1:nF

        fsw = f_sw_list(i);
        f0  = 2*fsw;

        % ---- Update design params for this frequency ----
        design_params_i = design_params;
        design_params_i.f = f0;   % adjust this field name if needed

        % ---- Run the three single-frequency functions ----
        out1 = analyzePriSwitches(Pmax, fsw, Ir_rms, 2, 10, ...
            [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 10000, ...
            'GaN Data tf.xlsx', 'SiC Data tf.xlsx');

        out2 = analyzeSecSwitches(Pmax, f0, 1, 10, [4 6 8], [4 6]);

        opt = optimize_VIRT(Pv_max_list, w_height_list, design_params_i);

        % ---- Efficiency calculation ----
        effOut = calcEfficiency_v3(out1, out2, "pareto", "pareto", ...
            Pmax, opt.P_total_min, opt.opt_design.A_footprint * 1e6, ...
            marker_by_jj, figNum + i, showEachTable, false);

        % ---- Save full result ----
        results(i).fsw   = fsw;
        results(i).out1  = out1;
        results(i).out2  = out2;
        results(i).opt   = opt;
        results(i).effOut = effOut;

        % ---- Use overall Pareto points if available ----
        if isfield(effOut, 'overallParetoPoints') && ~isempty(effOut.overallParetoPoints)
            ptsThis = effOut.overallParetoPoints;
        else
            ptsThis = effOut.points;
        end

        % ---- Append frequency tag to each point ----
        for k = 1:numel(ptsThis)
            allPts(end+1).fsw      = fsw; 
            allPts(end).area_mm2   = ptsThis(k).area_mm2;
            allPts(end).loss_W     = ptsThis(k).loss_W;
            allPts(end).eta        = ptsThis(k).eta;
            allPts(end).jj_pri     = ptsThis(k).jj_pri;
            allPts(end).jj_sec     = ptsThis(k).jj_sec;
            allPts(end).pri        = ptsThis(k).pri;
            allPts(end).sec        = ptsThis(k).sec;
        end

    end

    % ---------------- Build summary table ----------------
    nPts = numel(allPts);

    fsw_col    = reshape([allPts.fsw], [], 1);
    area_mm2   = reshape([allPts.area_mm2], [], 1);
    area_in2   = area_mm2 * 0.00155;
    loss_W     = reshape([allPts.loss_W], [], 1);
    eta        = reshape([allPts.eta], [], 1);
    jj_pri     = reshape([allPts.jj_pri], [], 1);
    jj_sec     = reshape([allPts.jj_sec], [], 1);

    pri_dev = strings(nPts,1);
    sec_dev = strings(nPts,1);
    for k = 1:nPts
        pri_dev(k) = string(allPts(k).pri.device_name);
        sec_dev(k) = string(allPts(k).sec.device_name);
    end

    point_id = (1:nPts)';

    T = table(point_id, fsw_col, eta, loss_W, area_mm2, area_in2, ...
        pri_dev, jj_pri, sec_dev, jj_sec, ...
        'VariableNames', {'point_id','fsw_Hz','Efficiency','Total_Loss_W', ...
        'Area_mm2','Area_in2','pri_device','Num_Parallel_Pri', ...
        'sec_device','Num_Parallel_Sec'});

    % ---------------- Combined comparison plot ----------------
    plotParetoAcrossFreq_v1(allPts, marker_by_jj, figNum, false);
    plotParetoAcrossFreq_v1(allPts, marker_by_jj, figNum + 1, true);

    % ---------------- Output ----------------
    freqSweepOut = struct();
    freqSweepOut.f_sw_list = f_sw_list;
    freqSweepOut.results = results;
    freqSweepOut.allParetoPoints = allPts;
    freqSweepOut.table = T;

end

function plotParetoAcrossFreq_v1(allPts, marker_by_jj, figNum, useIn2)

figure(figNum); clf;
hold on; grid on;

N = numel(allPts);
msize = 90;

% --------------------------------------------------
% Extract primary parallel counts and frequencies
% --------------------------------------------------
jj_pri = reshape([allPts.jj_pri], [], 1);
fswValsEach = reshape([allPts.fsw], [], 1);

jjPriVals = unique(jj_pri,'stable');
freqVals  = unique(fswValsEach,'stable');

% Color map for frequency
colors = lines(numel(freqVals));

% --------------------------------------------------
% Plot points
% --------------------------------------------------
for k = 1:N

    % ----- marker = primary parallel count -----
    jj_p = allPts(k).jj_pri;

    mk = 'o';
    if jj_p <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_p})
        mk = marker_by_jj{jj_p};
    end

    % ----- color = frequency -----
    fsw = allPts(k).fsw;
    colorIdx = find(freqVals == fsw,1);
    Color = colors(colorIdx,:);

    x = allPts(k).area_mm2;
    if useIn2
        x = x * 0.00155;
    end

    y = allPts(k).loss_W;

    scatter(x, y, msize, ...
        'Marker', mk, ...
        'MarkerFaceColor', Color, ...
        'MarkerEdgeColor', Color, ...
        'LineWidth', 1.2);

    text(x, y, sprintf('  %.2f%%',100*allPts(k).eta), ...
        'FontSize',10, ...
        'VerticalAlignment','middle');
end

% --------------------------------------------------
% Axis labels
% --------------------------------------------------
if useIn2
    xlabel('Total Area [in$^2$]','Interpreter','latex','FontSize',15);
else
    xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
end

ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);
title('Pareto Comparison Across Frequency','FontSize',15);

% --------------------------------------------------
% Marker legend (primary parallel count)
% --------------------------------------------------
h1 = gobjects(numel(jjPriVals),1);

for i = 1:numel(jjPriVals)

    jj = jjPriVals(i);

    mk = 'o';
    if jj <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj})
        mk = marker_by_jj{jj};
    end

    h1(i) = scatter(nan,nan,90, ...
        'Marker',mk, ...
        'MarkerFaceColor','k', ...
        'MarkerEdgeColor','k', ...
        'LineWidth',1.2);
end

% --------------------------------------------------
% Color legend (frequency)
% --------------------------------------------------
h2 = gobjects(numel(freqVals),1);

for i = 1:numel(freqVals)

    h2(i) = scatter(nan,nan,90, ...
        'Marker','o', ...
        'MarkerFaceColor',colors(i,:), ...
        'MarkerEdgeColor',colors(i,:), ...
        'LineWidth',1.2);
end

% Build frequency labels
freqLabels = strings(numel(freqVals),1);
for i = 1:numel(freqVals)
    freqLabels(i) = sprintf('f_{sw} = %.3g kHz', freqVals(i)/1e3);
end

legend([h1;h2], ...
       [compose('Pri %d',jjPriVals); freqLabels], ...
       'Location','bestoutside', ...
       'FontSize',10);

end