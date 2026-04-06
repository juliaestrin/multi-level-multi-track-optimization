function effOut = calcEfficiency_v4(PriOut, SecOut, modePri, modeSec, Power, other_loss, other_area, marker_by_jj, figNum, showTable, doPlot)
% calcEfficiency_v4
% Cases:
%   1) neither side is "pareto" -> scalar result
%   2) exactly one side is "pareto" -> point list + plots + table
%   3) both sides are "pareto" -> all pairwise combinations + plots + table
%
% Plot rule:
%   marker = primary-side jj25 for multilevel primary, otherwise primary jj
%   color  = secondary parallel count
%
% Point labels:
%   - efficiency shown to the right of the point
%   - jj75 shown above the point for multilevel primary only
%
% Primary-side point storage:
%   - normal topology:
%         jj_pri   = scalar jj
%         jj25_pri = NaN
%         jj75_pri = NaN
%   - multilevel topology:
%         jj_pri   = NaN
%         jj25_pri = value
%         jj75_pri = value

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
                 'jj_pri',{},'jj25_pri',{},'jj75_pri',{},'jj_sec',{}, ...
                 'pri',{},'sec',{});

    if isParetoPri && ~isParetoSec
        % ---------------- Case B1: pri pareto, sec fixed ----------------
        paretoSide = "pri";
        priCands = PriOut.best.pareto;
        secFixed = SecOut.best.(modeSec);

        N = numel(priCands);
        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[], ...
                            'jj_pri',[],'jj25_pri',[],'jj75_pri',[],'jj_sec',[], ...
                            'pri',[],'sec',[]), N, 1);

        for k = 1:N
            loss_total = priCands(k).loss_W + secFixed.loss_W + other_loss;
            area_total = priCands(k).area_mm2 + secFixed.area_mm2 + other_area;
            eta = (Power - loss_total) / Power;

            [jj_pri, jj25_pri, jj75_pri] = getPriParallelFields(priCands(k));
            jj_sec = getSecParallelField(secFixed);

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj_pri   = jj_pri;
            pts(k).jj25_pri = jj25_pri;
            pts(k).jj75_pri = jj75_pri;
            pts(k).jj_sec   = jj_sec;
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
                            'jj_pri',[],'jj25_pri',[],'jj75_pri',[],'jj_sec',[], ...
                            'pri',[],'sec',[]), N, 1);

        for k = 1:N
            loss_total = priFixed.loss_W + secCands(k).loss_W + other_loss;
            area_total = priFixed.area_mm2 + secCands(k).area_mm2 + other_area;
            eta = (Power - loss_total) / Power;

            [jj_pri, jj25_pri, jj75_pri] = getPriParallelFields(priFixed);
            jj_sec = getSecParallelField(secCands(k));

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj_pri   = jj_pri;
            pts(k).jj25_pri = jj25_pri;
            pts(k).jj75_pri = jj75_pri;
            pts(k).jj_sec   = jj_sec;
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
                            'jj_pri',[],'jj25_pri',[],'jj75_pri',[],'jj_sec',[], ...
                            'pri',[],'sec',[]), N, 1);

        idx = 0;
        for i = 1:Np
            for j = 1:Ns
                idx = idx + 1;

                loss_total = priCands(i).loss_W + secCands(j).loss_W + other_loss;
                area_total = priCands(i).area_mm2 + secCands(j).area_mm2 + other_area;
                eta = (Power - loss_total) / Power;

                [jj_pri, jj25_pri, jj75_pri] = getPriParallelFields(priCands(i));
                jj_sec = getSecParallelField(secCands(j));

                pts(idx).loss_W   = loss_total;
                pts(idx).area_mm2 = area_total;
                pts(idx).eta      = eta;
                pts(idx).jj_pri   = jj_pri;
                pts(idx).jj25_pri = jj25_pri;
                pts(idx).jj75_pri = jj75_pri;
                pts(idx).jj_sec   = jj_sec;
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

    pri_parallel = strings(nPts,1);
    sec_parallel = strings(nPts,1);

    for k = 1:nPts
        pri_parallel(k) = makePriParallelString(pts(k));
        sec_parallel(k) = sprintf('%d', pts(k).jj_sec);
    end

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
              pri_dev, pri_parallel, pri_loss_W, pri_area_mm2, ...
              sec_dev, sec_parallel, sec_loss_W, sec_area_mm2, ...
              other_loss_col, other_area_col, ...
              pri_area_mm2./area_mm2, sec_area_mm2./area_mm2, other_area_col./area_mm2, ...
              pri_loss_W./loss_W, sec_loss_W./loss_W, other_loss_col./loss_W, ...
              'VariableNames', {'point_id','Efficiency','Total Loss [W]','Area [mm2]','Area [in2]', ...
                                'pri_device','Primary Parallel','pri_loss_W','pri_area_mm2', ...
                                'sec_device','Secondary Parallel','sec_loss_W','sec_area_mm2', ...
                                'other_loss_W','other_area_mm2','% area pri','% area sec','% area mag', ...
                                '% loss pri','% loss sec','% loss mag'});

    T = sortrows(T, {'Total Loss [W]','Area [mm2]'}, {'ascend','ascend'});
    effOut.table = T;

    if showTable
        disp(T);

        if doPlot
            fTbl = figure(figNum + 2); clf(fTbl);
            fTbl.Name = sprintf('Pareto Results Table (%s pareto)', paretoSide);
            uitable(fTbl, 'Data', T{:,:}, 'ColumnName', T.Properties.VariableNames, ...
                    'Units','normalized', 'Position',[0 0 1 1]);
        end
    end

    % ---------------- Overall Pareto front data for pareto/pareto case ----------------
    if isParetoPri && isParetoSec && numel(SecOut.best.pareto) > 1

        area_vec = reshape([pts.area_mm2], [], 1);
        loss_vec = reshape([pts.loss_W],   [], 1);

        isPareto = computeParetoFront(area_vec, loss_vec);
        paretoIdx = find(isPareto);

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
        plotAreaLoss_v4(pts, marker_by_jj, figNum,   false);
        plotAreaLoss_v4(pts, marker_by_jj, figNum+1, true);

        if ~isempty(effOut.overallParetoPoints)
            plotOverallParetoFront_v4(pts, effOut.overallParetoIdx, marker_by_jj, figNum+3, false);
            plotOverallParetoFront_v4(pts, effOut.overallParetoIdx, marker_by_jj, figNum+4, true);
        end
    end

end


function [jj_pri, jj25_pri, jj75_pri] = getPriParallelFields(one)
    if isfield(one, 'jj')
        jj_pri   = one.jj;
        jj25_pri = NaN;
        jj75_pri = NaN;
    elseif isfield(one, 'jj25') && isfield(one, 'jj75')
        jj_pri   = NaN;
        jj25_pri = one.jj25;
        jj75_pri = one.jj75;
    else
        error('Primary candidate does not contain jj or (jj25,jj75).');
    end
end


function jj_sec = getSecParallelField(one)
    if isfield(one, 'jj')
        jj_sec = one.jj;
    else
        error('Secondary candidate does not contain scalar jj.');
    end
end


function s = makePriParallelString(pt)
    if ~isnan(pt.jj_pri)
        s = sprintf('%d', pt.jj_pri);
    else
        s = sprintf('(%d,%d)', pt.jj25_pri, pt.jj75_pri);
    end
end


function mk = getPrimaryMarker(pt, marker_by_jj)
    if ~isnan(pt.jj_pri)
        markerKey = pt.jj_pri;
    else
        markerKey = pt.jj25_pri;   % requested rule
    end

    mk = 'o';
    if markerKey <= numel(marker_by_jj) && ~isempty(marker_by_jj{markerKey})
        mk = marker_by_jj{markerKey};
    end
end


function plotAreaLoss_v4(pts, marker_by_jj, figNum, useIn2)

    figure(figNum); clf;
    hold on; grid on;

    N = numel(pts);
    msize = 90;

    % --------------------------------------------------
    % Extract primary / secondary info
    % --------------------------------------------------
    secVals = reshape([pts.jj_sec], [], 1);
    secUnique = unique(secVals, 'stable');
    colors = lines(numel(secUnique));

    priLegendTexts = strings(N,1);
    priMarkerKeys  = nan(N,1);

    for k = 1:N
        priLegendTexts(k) = makePriParallelString(pts(k));

        if ~isnan(pts(k).jj_pri)
            priMarkerKeys(k) = pts(k).jj_pri;
        else
            priMarkerKeys(k) = pts(k).jj25_pri;
        end
    end

    [priLegendUnique, iaPri] = unique(priLegendTexts, 'stable');
    priMarkerUnique = priMarkerKeys(iaPri);

    % --------------------------------------------------
    % Plot points
    % --------------------------------------------------
    for k = 1:N
        mk = getPrimaryMarker(pts(k), marker_by_jj);

        secIdx = find(secUnique == pts(k).jj_sec, 1);
        Color = colors(secIdx,:);

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

        % Efficiency on the right
        text(x, y, sprintf('  %.2f%%', 100*pts(k).eta), ...
            'FontSize',10, ...
            'VerticalAlignment','middle', ...
            'HorizontalAlignment','left');

        % jj75 above for multilevel-primary points only
        % if isnan(pts(k).jj_pri)
        %     text(x, y, sprintf('%d', pts(k).jj75_pri), ...
        %         'FontSize',10, ...
        %         'VerticalAlignment','bottom', ...
        %         'HorizontalAlignment','center');
        % end
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
    % Marker legend (primary topology setting)
    % --------------------------------------------------
    h1 = gobjects(numel(priLegendUnique),1);

    for i = 1:numel(priLegendUnique)
        mk = 'o';
        key = priMarkerUnique(i);
        if key <= numel(marker_by_jj) && ~isempty(marker_by_jj{key})
            mk = marker_by_jj{key};
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
    h2 = gobjects(numel(secUnique),1);

    for i = 1:numel(secUnique)
        h2(i) = scatter(nan,nan,90, ...
            'Marker','o', ...
            'MarkerFaceColor',colors(i,:), ...
            'MarkerEdgeColor',colors(i,:), ...
            'LineWidth',1.2);
    end

    legend([h1;h2], ...
           ["Pri " + priLegendUnique; compose('Sec %d',secUnique)], ...
           'Location','bestoutside', ...
           'FontSize',10);

end


function plotOverallParetoFront_v4(pts, paretoIdx, marker_by_jj, figNum, useIn2)

    pfPts = pts(paretoIdx);

    figure(figNum); clf;
    hold on; grid on;

    N = numel(pfPts);
    msize = 90;

    secVals = reshape([pfPts.jj_sec], [], 1);
    secUnique = unique(secVals, 'stable');
    colors = lines(numel(secUnique));

    priLegendTexts = strings(N,1);
    priMarkerKeys  = nan(N,1);

    for k = 1:N
        priLegendTexts(k) = makePriParallelString(pfPts(k));

        if ~isnan(pfPts(k).jj_pri)
            priMarkerKeys(k) = pfPts(k).jj_pri;
        else
            priMarkerKeys(k) = pfPts(k).jj25_pri;
        end
    end

    [priLegendUnique, iaPri] = unique(priLegendTexts, 'stable');
    priMarkerUnique = priMarkerKeys(iaPri);

    % --------------------------------------------------
    % Plot Pareto-front points only
    % --------------------------------------------------
    for k = 1:N
        mk = getPrimaryMarker(pfPts(k), marker_by_jj);

        secIdx = find(secUnique == pfPts(k).jj_sec, 1);
        Color = colors(secIdx,:);

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

        % Efficiency on the right
        text(x, y, sprintf('  %.2f%%', 100*pfPts(k).eta), ...
            'FontSize',10, ...
            'VerticalAlignment','middle', ...
            'HorizontalAlignment','left');

        % jj75 above for multilevel-primary points only
        if isnan(pfPts(k).jj_pri)
            text(x, y, sprintf('j_{75}=%d', pfPts(k).jj75_pri), ...
                'FontSize',10, ...
                'VerticalAlignment','bottom', ...
                'HorizontalAlignment','center');
        end
    end

    if useIn2
        xlabel('Total Area [in$^2$]','Interpreter','latex','FontSize',15);
    else
        xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
    end

    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);
    title('Overall Pareto Front','FontSize',15);

    % --------------------------------------------------
    % Marker legend (primary topology setting)
    % --------------------------------------------------
    h1 = gobjects(numel(priLegendUnique),1);

    for i = 1:numel(priLegendUnique)
        mk = 'o';
        key = priMarkerUnique(i);
        if key <= numel(marker_by_jj) && ~isempty(marker_by_jj{key})
            mk = marker_by_jj{key};
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
    h2 = gobjects(numel(secUnique),1);

    for i = 1:numel(secUnique)
        h2(i) = scatter(nan,nan,90, ...
            'Marker','o', ...
            'MarkerFaceColor',colors(i,:), ...
            'MarkerEdgeColor',colors(i,:), ...
            'LineWidth',1.2);
    end

    legend([h1;h2], ...
           ["Pri " + priLegendUnique; compose('Sec %d',secUnique)], ...
           'Location','bestoutside', ...
           'FontSize',10);

end