function effOut = calcEfficiency_v3(PriOut, SecOut, modePri, modeSec, Power, other_loss, other_area, marker_by_jj, figNum, showTable, doPlot)
% calcEfficiency_v3
% Cases:
%   1) neither side is "pareto" -> scalar result
%   2) exactly one side is "pareto" -> point list + plots + table
%   3) both sides are "pareto" -> all pairwise combinations + plots + table
%
% Plot rule:
%   marker = primary parallel count
%   color  = secondary parallel count

    % ---------------- Defaults ----------------
    if nargin < 6 || isempty(other_loss), other_loss = 0; end
    if nargin < 7 || isempty(other_area), other_area = 0; end
    if nargin < 8 || isempty(marker_by_jj), marker_by_jj = makeMarkerMap(); end
    if nargin < 9 || isempty(figNum), figNum = 10; end
    if nargin < 10 || isempty(showTable), showTable = true; end
    if nargin < 11 || isempty(doPlot), doPlot = true; end

    if ~isstring(modePri) || ~isscalar(modePri) || ~isstring(modeSec) || ~isscalar(modeSec)
        error('modePri and modeSec must be string scalars (e.g., "minLoss", "pareto").');
    end

    isParetoPri = (modePri == "pareto");
    isParetoSec = (modeSec == "pareto");

    % ---------------- Case A: no pareto -> scalar output ----------------
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
        effOut.table = table();
        return;
    end

    % ============================================================
    % Build points for pareto cases
    % ============================================================
    pts = struct('area_mm2',{},'loss_W',{},'eta',{}, ...
                 'jj_pri',{},'jj_sec',{}, ...
                 'pri',{},'sec',{});

    if isParetoPri && ~isParetoSec
        % ---------------- Case B1: pri pareto, sec fixed ----------------
        paretoSide = "pri";
        priCands = PriOut.best.pareto;
        secFixed = SecOut.best.(modeSec);

        N = numel(priCands);
        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[], ...
                            'jj_pri',[],'jj_sec',[], ...
                            'pri',[],'sec',[]), N, 1);

        for k = 1:N
            loss_total = priCands(k).loss_W + secFixed.loss_W + other_loss;
            area_total = priCands(k).area_mm2 + secFixed.area_mm2 + other_area;
            eta = (Power - loss_total) / Power;

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj_pri   = priCands(k).jj;
            pts(k).jj_sec   = secFixed.jj;
            pts(k).pri      = priCands(k);
            pts(k).sec      = secFixed;
        end

    elseif ~isParetoPri && isParetoSec
        % ---------------- Case B2: sec pareto, pri fixed ----------------
        paretoSide = "sec";
        secCands = SecOut.best.pareto;
        priFixed = PriOut.best.(modePri);

        N = numel(secCands);
        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[], ...
                            'jj_pri',[],'jj_sec',[], ...
                            'pri',[],'sec',[]), N, 1);

        for k = 1:N
            loss_total = priFixed.loss_W + secCands(k).loss_W + other_loss;
            area_total = priFixed.area_mm2 + secCands(k).area_mm2 + other_area;
            eta = (Power - loss_total) / Power;

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj_pri   = priFixed.jj;
            pts(k).jj_sec   = secCands(k).jj;
            pts(k).pri      = priFixed;
            pts(k).sec      = secCands(k);
        end

    else
        % ---------------- Case C: both pri and sec are pareto ----------------
        paretoSide = "both";
        priCands = PriOut.best.pareto;
        secCands = SecOut.best.pareto;

        Np = numel(priCands);
        Ns = numel(secCands);
        N  = Np * Ns;

        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[], ...
                            'jj_pri',[],'jj_sec',[], ...
                            'pri',[],'sec',[]), N, 1);

        idx = 0;
        for i = 1:Np
            for j = 1:Ns
                idx = idx + 1;

                loss_total = priCands(i).loss_W + secCands(j).loss_W + other_loss;
                area_total = priCands(i).area_mm2 + secCands(j).area_mm2 + other_area;
                eta = (Power - loss_total) / Power;

                pts(idx).loss_W   = loss_total;
                pts(idx).area_mm2 = area_total;
                pts(idx).eta      = eta;
                pts(idx).jj_pri   = priCands(i).jj;
                pts(idx).jj_sec   = secCands(j).jj;
                pts(idx).pri      = priCands(i);
                pts(idx).sec      = secCands(j);
            end
        end
    end

    effOut = struct();
    effOut.points = pts;
    effOut.paretoSide = paretoSide;

    % ---------------- Build RESULTS TABLE ----------------
    nPts = numel(pts);

    area_mm2 = reshape([pts.area_mm2], [], 1);
    area_in2 = area_mm2 * 0.00155;
    loss_W   = reshape([pts.loss_W],   [], 1);
    eta      = reshape([pts.eta],      [], 1);

    pri_loss_W   = reshape(arrayfun(@(p) p.pri.loss_W,   pts), [], 1);
    sec_loss_W   = reshape(arrayfun(@(p) p.sec.loss_W,   pts), [], 1);
    pri_area_mm2 = reshape(arrayfun(@(p) p.pri.area_mm2, pts), [], 1);
    sec_area_mm2 = reshape(arrayfun(@(p) p.sec.area_mm2, pts), [], 1);

    jj_pri = reshape([pts.jj_pri], [], 1);
    jj_sec = reshape([pts.jj_sec], [], 1);

    pri_dev = strings(nPts,1);
    sec_dev = strings(nPts,1);
    for k = 1:nPts
        pri_dev(k) = string(pts(k).pri.device_name);
        sec_dev(k) = string(pts(k).sec.device_name);
    end

    point_id = (1:nPts)';

    other_loss_col = repmat(other_loss, nPts, 1);
    other_area_col = repmat(other_area, nPts, 1);

    T = table(point_id, eta, loss_W, area_mm2, area_in2, ...
              pri_dev, jj_pri, pri_loss_W, pri_area_mm2, ...
              sec_dev, jj_sec, sec_loss_W, sec_area_mm2, ...
              other_loss_col, other_area_col, ...
              pri_area_mm2./area_mm2, sec_area_mm2./area_mm2, other_area_col./area_mm2, ...
              pri_loss_W./loss_W, sec_loss_W./loss_W, other_loss_col./loss_W, ...
              'VariableNames', {'point_id','Efficiency','Total Loss [W]','Area [mm2]','Area [in2]', ...
                                'pri_device','Num Parallel Pri','pri_loss_W','pri_area_mm2', ...
                                'sec_device','Num Parallel Sec','sec_loss_W','sec_area_mm2', ...
                                'other_loss_W','other_area_mm2','% area pri','% area sec','% area mag', ...
                                '% loss pri','% loss sec','% loss mag'});

    T = sortrows(T, {'Total Loss [W]','Area [mm2]'}, {'ascend','ascend'});
    effOut.table = T;

    if showTable
        disp(T);

        if doPlot
            fTbl = figure(figNum + 2); clf(fTbl);
            fTbl.Name = sprintf('Pareto Results Table (%s pareto)', paretoSide);
            uitable(fTbl, 'Data', T, 'Units','normalized', 'Position',[0 0 1 1]);
        end
    end

    % ---------------- Overall Pareto front data for pareto/pareto case ----------------
    % Only generate/store this when the secondary Pareto set has > 1 points
    if isParetoPri && isParetoSec && numel(SecOut.best.pareto) > 1

        area_vec = reshape([pts.area_mm2], [], 1);
        loss_vec = reshape([pts.loss_W],   [], 1);

        isPareto = computeParetoFront(area_vec, loss_vec);
        paretoIdx = find(isPareto);

        % Sort by area for cleaner visual ordering
        [~, order] = sort(area_vec(paretoIdx), 'ascend');
        paretoIdx = paretoIdx(order);

        effOut.overallParetoIdx = paretoIdx;
        effOut.overallParetoPoints = pts(paretoIdx);

    else
        effOut.overallParetoIdx = [];
        effOut.overallParetoPoints = struct([]);
    end

    % ---------------- Plotting ----------------
    if doPlot
        plotAreaLoss_v3(pts, marker_by_jj, figNum,   false);
        plotAreaLoss_v3(pts, marker_by_jj, figNum+1, true);

        if ~isempty(effOut.overallParetoPoints)
            plotOverallParetoFront_v3(pts, effOut.overallParetoIdx, marker_by_jj, figNum+3, false);
            plotOverallParetoFront_v3(pts, effOut.overallParetoIdx, marker_by_jj, figNum+4, true);
        end
    end

end


function plotAreaLoss_v3(pts, marker_by_jj, figNum, useIn2)

figure(figNum); clf;
hold on; grid on;

N = numel(pts);
msize = 90;

% --------------------------------------------------
% Extract primary/secondary info
% --------------------------------------------------
jj_pri = reshape([pts.jj_pri], [], 1);

sec_names = strings(N,1);
for k = 1:N
    sec_names(k) = string(pts(k).sec.device_name);
end

% Unique sets
jjPriVals = unique(jj_pri,'stable');
secNameVals = unique(sec_names,'stable');

% Colors for secondary points
colors = lines(numel(secNameVals));

% --------------------------------------------------
% Plot points
% --------------------------------------------------
for k = 1:N

    % ----- marker = primary parallel -----
    jj_p = pts(k).jj_pri;

    mk = 'o';
    if jj_p <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_p})
        mk = marker_by_jj{jj_p};
    end

    % ----- color = secondary candidate -----
    secName = string(pts(k).sec.device_name);
    colorIdx = find(secNameVals == secName,1);
    Color = colors(colorIdx,:);

    x = pts(k).area_mm2;
    if useIn2
        x = x * 0.00155;
    end

    y = pts(k).loss_W;

    scatter(x, y, msize, ...
        'Marker', mk, ...
        'MarkerFaceColor', Color, ...
        'MarkerEdgeColor', Color, ...
        'LineWidth', 1.2);

    text(x, y, sprintf('  %.2f%%',100*pts(k).eta), ...
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

% --------------------------------------------------
% Marker legend (primary parallel)
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
% Color legend (secondary candidate)
% --------------------------------------------------
h2 = gobjects(numel(secNameVals),1);

for i = 1:numel(secNameVals)

    h2(i) = scatter(nan,nan,90, ...
        'Marker','o', ...
        'MarkerFaceColor',colors(i,:), ...
        'MarkerEdgeColor',colors(i,:), ...
        'LineWidth',1.2);
end

legend([h1;h2], ...
       [compose('Pri %d',jjPriVals); secNameVals], ...
       'Location','bestoutside', ...
       'FontSize',10);

end

function plotOverallParetoFront_v3(pts, paretoIdx, marker_by_jj, figNum, useIn2)

pfPts = pts(paretoIdx);

figure(figNum); clf;
hold on; grid on;

N = numel(pfPts);
msize = 90;

% --------------------------------------------------
% Extract primary/secondary parallel counts
% --------------------------------------------------
jj_pri = reshape([pfPts.jj_pri], [], 1);
jj_sec = reshape([pfPts.jj_sec], [], 1);

% Unique sets
jjPriVals = unique(jj_pri,'stable');
jjSecVals = unique(jj_sec,'stable');

% Colors for secondary parallel counts
colors = lines(numel(jjSecVals));

% --------------------------------------------------
% Plot Pareto-front points only
% --------------------------------------------------
for k = 1:N

    % ----- marker = primary parallel count -----
    jj_p = pfPts(k).jj_pri;

    mk = 'o';
    if jj_p <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_p})
        mk = marker_by_jj{jj_p};
    end

    % ----- color = secondary parallel count -----
    jj_s = pfPts(k).jj_sec;
    colorIdx = find(jjSecVals == jj_s, 1);
    Color = colors(colorIdx,:);

    x = pfPts(k).area_mm2;
    if useIn2
        x = x * 0.00155;
    end

    y = pfPts(k).loss_W;

    scatter(x, y, msize, ...
        'Marker', mk, ...
        'MarkerFaceColor', Color, ...
        'MarkerEdgeColor', Color, ...
        'LineWidth', 1.2);

    text(x, y, sprintf('  %.2f%%',100*pfPts(k).eta), ...
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
title('Overall Pareto Front','FontSize',15);

% --------------------------------------------------
% Marker legend (primary parallel)
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
% Color legend (secondary parallel)
% --------------------------------------------------
h2 = gobjects(numel(jjSecVals),1);

for i = 1:numel(jjSecVals)
    h2(i) = scatter(nan,nan,90, ...
        'Marker','o', ...
        'MarkerFaceColor',colors(i,:), ...
        'MarkerEdgeColor',colors(i,:), ...
        'LineWidth',1.2);
end

legend([h1;h2], ...
       [compose('Pri %d',jjPriVals); compose('Sec %d',jjSecVals)], ...
       'Location','bestoutside', ...
       'FontSize',10);

end

function plotOverallParetoFront_old(pts, paretoIdx, marker_by_jj, figNum, useIn2)

figure(figNum); clf;
hold on; grid on;

N = numel(pts);
msizeAll = 70;
msizePf  = 120;

% --------------------------------------------------
% Extract primary / secondary parallel counts
% --------------------------------------------------
jj_pri = reshape([pts.jj_pri], [], 1);
jj_sec = reshape([pts.jj_sec], [], 1);

jjPriVals = unique(jj_pri, 'stable');
jjSecVals = unique(jj_sec, 'stable');

% Color map for secondary-side parallel count
colors = lines(numel(jjSecVals));

% --------------------------------------------------
% Plot all points in background
% --------------------------------------------------
for k = 1:N

    % marker = primary-side parallel count
    jj_p = pts(k).jj_pri;
    mk = 'o';
    if jj_p <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_p})
        mk = marker_by_jj{jj_p};
    end

    % color = secondary-side parallel count
    jj_s = pts(k).jj_sec;
    colorIdx = find(jjSecVals == jj_s, 1);
    Color = colors(colorIdx,:);

    x = pts(k).area_mm2;
    if useIn2
        x = x * 0.00155;
    end
    y = pts(k).loss_W;

    scatter(x, y, msizeAll, ...
        'Marker', mk, ...
        'MarkerFaceColor', 'none', ...
        'MarkerEdgeColor', Color, ...
        'LineWidth', 0.8);
end

% --------------------------------------------------
% Pareto front coordinates
% --------------------------------------------------
pfPts = pts(paretoIdx);

xpf = reshape([pfPts.area_mm2], [], 1);
if useIn2
    xpf = xpf * 0.00155;
end
ypf = reshape([pfPts.loss_W], [], 1);

plot(xpf, ypf, '-', 'LineWidth', 1.8, 'Color', 'k');

% --------------------------------------------------
% Highlight Pareto-front points
% --------------------------------------------------
for ii = 1:numel(paretoIdx)

    k = paretoIdx(ii);

    % marker = primary-side parallel count
    jj_p = pts(k).jj_pri;
    mk = 'o';
    if jj_p <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_p})
        mk = marker_by_jj{jj_p};
    end

    % color = secondary-side parallel count
    jj_s = pts(k).jj_sec;
    colorIdx = find(jjSecVals == jj_s, 1);
    Color = colors(colorIdx,:);

    x = pts(k).area_mm2;
    if useIn2
        x = x * 0.00155;
    end
    y = pts(k).loss_W;

    scatter(x, y, msizePf, ...
        'Marker', mk, ...
        'MarkerFaceColor', Color, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.4);

    text(x, y, sprintf('  %.2f%%', 100*pts(k).eta), ...
        'FontSize', 10, ...
        'VerticalAlignment', 'middle');
end

% --------------------------------------------------
% Labels
% --------------------------------------------------
if useIn2
    xlabel('Total Area [in$^2$]', 'Interpreter','latex', 'FontSize',15);
else
    xlabel('Total Area [mm$^2$]', 'Interpreter','latex', 'FontSize',15);
end
ylabel('Total Power Loss [W]', 'Interpreter','latex', 'FontSize',15);
title('Overall Pareto Front', 'FontSize', 15);

% --------------------------------------------------
% Legend
% --------------------------------------------------
h1 = gobjects(numel(jjPriVals),1);
for i = 1:numel(jjPriVals)
    jj = jjPriVals(i);

    mk = 'o';
    if jj <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj})
        mk = marker_by_jj{jj};
    end

    h1(i) = scatter(nan, nan, 90, ...
        'Marker', mk, ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.2);
end

h2 = gobjects(numel(jjSecVals),1);
for i = 1:numel(jjSecVals)
    h2(i) = scatter(nan, nan, 90, ...
        'Marker', 'o', ...
        'MarkerFaceColor', colors(i,:), ...
        'MarkerEdgeColor', colors(i,:), ...
        'LineWidth', 1.2);
end

h3 = plot(nan, nan, '-', 'LineWidth', 1.8, 'Color', 'k');

legend([h1; h2; h3], ...
       [compose('Pri %d', jjPriVals); compose('Sec %d', jjSecVals); "Overall Pareto front"], ...
       'Location', 'bestoutside', ...
       'FontSize', 10);

end