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