function effOut = calcEfficiency_v2(PriOut, SecOut, modePri, modeSec, Power, other_loss, other_area, marker_by_jj, figNum, showTable)
% calcEfficiency
% - If neither side is "pareto": returns scalar results (eta, total loss, total area)
% - If exactly one side is "pareto": builds list of points, plots, and creates a results table

    % ---------------- Defaults ----------------
    if nargin < 6 || isempty(other_loss), other_loss = 0; end
    if nargin < 7 || isempty(other_area), other_area = 0; end
    if nargin < 8 || isempty(marker_by_jj), marker_by_jj = makeMarkerMap(); end
    if nargin < 9 || isempty(figNum), figNum = 10; end
    if nargin < 10 || isempty(showTable), showTable = true; end

    if ~isstring(modePri) || ~isscalar(modePri) || ~isstring(modeSec) || ~isscalar(modeSec)
        error('modePri and modeSec must be string scalars (e.g., "minLoss", "pareto").');
    end

    isParetoPri = (modePri == "pareto");
    isParetoSec = (modeSec == "pareto");

    if isParetoPri && isParetoSec
        error('Both modePri and modeSec cannot be "pareto" in this function.');
    end

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
        effOut.table = table(); % empty
        return;
    end

    % ---------------- Case B: exactly one pareto -> compute list ----------------
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
            pts(k).jj       = priCands(k).jj;
            pts(k).pri      = priCands(k);
            pts(k).sec      = secFixed;
        end
    else
        N = numel(secCands);
        pts = repmat(struct('area_mm2',[],'loss_W',[],'eta',[],'jj',[],'pri',[],'sec',[]), N, 1);
        for k = 1:N
            loss_total = priFixed.loss_W + secCands(k).loss_W + other_loss;
            area_total = priFixed.area_mm2 + secCands(k).area_mm2 + other_area;  % <- include other_area (bugfix)
            eta = (Power - loss_total) / Power;

            pts(k).loss_W   = loss_total;
            pts(k).area_mm2 = area_total;
            pts(k).eta      = eta;
            pts(k).jj       = secCands(k).jj;
            pts(k).pri      = priFixed;
            pts(k).sec      = secCands(k);
        end
    end

    effOut = struct();
    effOut.points = pts;
    effOut.paretoSide = paretoSide;

    % ---------------- Build RESULTS TABLE ----------------
    % Core totals
    area_mm2 = [pts.area_mm2]';
    area_in2 = area_mm2 * 0.00155;
    loss_W   = [pts.loss_W]';
    eta      = [pts.eta]';
    eta_pct  = 100*eta;
    jj       = [pts.jj]';

    % Component details (totals and fixed/candidate attributes)
    pri_loss_W   = arrayfun(@(p) p.pri.loss_W, pts)';
    sec_loss_W   = arrayfun(@(p) p.sec.loss_W, pts)';
    pri_area_mm2 = arrayfun(@(p) p.pri.area_mm2, pts)';
    sec_area_mm2 = arrayfun(@(p) p.sec.area_mm2, pts)';

    pri_dev = strings(numel(pts),1);
    sec_dev = strings(numel(pts),1);
    sec_parallel = strings(numel(pts),1);
    for k = 1:numel(pts)
        pri_dev(k) = string(pts(k).pri.device_name);
        sec_dev(k) = string(pts(k).sec.device_name);
        sec_parallel(k) = pts(k).sec.jj;
    end



    % Optional: add an easy ID + sort order
    point_id = (1:numel(pts))';

    T = table(point_id, eta, loss_W, area_mm2, area_in2, ...
              pri_dev, jj, pri_loss_W.', pri_area_mm2.',...
              sec_dev, sec_parallel, sec_loss_W.', sec_area_mm2.', ...
              repmat(other_loss, numel(pts), 1), repmat(other_area, numel(pts), 1), ...
              pri_area_mm2'./area_mm2, sec_area_mm2'./area_mm2, repmat(other_area, numel(pts), 1)./area_mm2,...
              pri_loss_W'./loss_W, sec_loss_W'./loss_W, repmat(other_loss, numel(pts), 1)./loss_W,...
              'VariableNames', {'point_id','Efficiency','Total Loss [W]','Area [mm2]','Area [in2]', ...
                                'pri_device', 'Num Parallel Pri', 'pri_loss_W', 'pri_area_mm2', ...
                                'sec_device', 'Num Parallel Sec','sec_loss_W','sec_area_mm2', ...
                                'other_loss_W','other_area_mm2','% area pri','% area sec','% area mag' ...
                                '% loss pri', '% loss sec','% loss mag'});

    % Sort (pick what you like). This is nice for reading:
    % primary: loss ascending, secondary: area ascending
    T = sortrows(T, {'Total Loss [W]','Area [mm2]'}, {'ascend','ascend'});

    effOut.table = T;

    if showTable
        disp(T);
        % Also pop a UI table (comment out if you don't want a window)
        fTbl = figure(); clf(fTbl);
        fTbl.Name = sprintf('Pareto Results Table');
        uitable(fTbl, 'Data', T, 'Units','normalized', 'Position',[0 0 1 1]);
    end

    % ---------------- Plotting (mm^2 and in^2) ----------------
    plotAreaLoss(pts, paretoSide, marker_by_jj, figNum,   false); % mm^2
    plotAreaLoss(pts, paretoSide, marker_by_jj, figNum+1, true);  % in^2
end


function plotAreaLoss(pts, paretoSide, marker_by_jj, figNum, useIn2)
    figure(); clf;
    hold on; grid on;

    N = numel(pts);
    msize = 90;

    % color by device (on pareto side)
    devNames = strings(N,1);
    for k = 1:N
        if paretoSide == "pri"
            devNames(k) = string(pts(k).pri.device_name);
        else
            devNames(k) = string(pts(k).sec.device_name);
        end
    end
    [~, ~, devIdx] = unique(devNames, 'stable');
    colors = lines(numel(unique(devNames,'stable')));

    for k = 1:N
        jj = pts(k).jj;

        mk = 'o';
        if jj <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj})
            mk = marker_by_jj{jj};
        end

        Color = colors(devIdx(k),:);

        x = pts(k).area_mm2;
        if useIn2, x = x * 0.00155; end
        y = pts(k).loss_W;

        scatter(x, y, msize, 'Marker', mk, ...
            'MarkerFaceColor', Color, 'MarkerEdgeColor', Color, 'LineWidth', 1.2);

        text(x, y, sprintf('  %.2f%%', 100*pts(k).eta), ...
            'FontSize', 10, 'VerticalAlignment', 'middle');
    end

    if useIn2
        xlabel('Total Area [in$^2$]','Interpreter','latex','FontSize',15);
    else
        xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
    end
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);

    jjVals = unique([pts.jj]);
    h = gobjects(numel(jjVals),1);
    for i = 1:numel(jjVals)
        jj = jjVals(i);
        mk = 'o';
        if jj <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj})
            mk = marker_by_jj{jj};
        end
        h(i) = scatter(nan, nan, 90, 'Marker', mk, ...
            'MarkerFaceColor', 'k', 'MarkerEdgeColor','k', 'LineWidth', 1.2);
    end
    legend(h, compose('%d-parallel', jjVals), 'Location','best', ...
        'Orientation','horizontal','FontSize',10);
end