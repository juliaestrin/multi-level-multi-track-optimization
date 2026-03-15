% Comparing to v1: here is the one with loss break down: P_cond, P_off and
% P_gate 
function best = pickBestDevices_v2(AreaFull, LossFull, LossAreaFull, PcondFull, PoffFull, PgateFull, pri_sw_count, jj_set, compare_list, Meta)
    cols = zeros(1, numel(compare_list));
    for k = 1:numel(compare_list)
        cols(k) = find(jj_set == compare_list(k), 1);
    end

    Loss = LossFull(:, cols);
    Area = AreaFull(:, cols);
    LossArea = LossAreaFull(:, cols);

    % Scale to total primary-side loss so they are consistent with Loss = P_total_plot
    Pcond = pri_sw_count * PcondFull(:, cols);
    Poff  = pri_sw_count * PoffFull(:, cols);
    Pgate = pri_sw_count * PgateFull(:, cols);

    feas = ~isnan(Loss);

    % Optional but recommended: keep only feasible points
    Pcond(~feas) = NaN;
    Poff(~feas)  = NaN;
    Pgate(~feas) = NaN;

    Loss1 = Loss;
    Loss1(~feas) = NaN;
    [~, lin] = min(Loss1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Loss1), lin);
    jj_abs = compare_list(colLocal);
    best.minLoss = packOneLocal_v2(dev, jj_abs, compare_list, Loss, Area, Pcond, Poff, Pgate, Meta);

    Area1 = Area;
    Area1(~feas) = NaN;
    [~, lin] = min(Area1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Area1), lin);
    jj_abs = compare_list(colLocal);
    best.minArea = packOneLocal_v2(dev, jj_abs, compare_list, Loss, Area, Pcond, Poff, Pgate, Meta);

    LA1 = LossArea;
    LA1(~feas) = NaN;
    [~, lin] = min(LA1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(LA1), lin);
    jj_abs = compare_list(colLocal);
    best.minLossArea = packOneLocal_v2(dev, jj_abs, compare_list, Loss, Area, Pcond, Poff, Pgate, Meta);

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
        'loss_W', [], 'area_mm2', [], 'loss_area_W_mm2_product', [], ...
        'P_cond_W', [], 'P_off_W', [], 'P_gate_W', []), ...
        numel(pDevIdx), 1);

    for m = 1:numel(pDevIdx)
        di = pDevIdx(m);
        cj = pColIdx(m);

        jj_abs = compare_list(cj);

        paretoList(m).device_index = di;
        paretoList(m).device_name  = string(Meta.Name(di));
        paretoList(m).jj           = jj_abs;

        paretoList(m).loss_W  = Loss(di, cj);
        paretoList(m).area_mm2 = Area(di, cj);
        paretoList(m).loss_area_W_mm2_product = Area(di, cj) * Loss(di, cj);

        paretoList(m).P_cond_W = Pcond(di, cj);
        paretoList(m).P_off_W  = Poff(di, cj);
        paretoList(m).P_gate_W = Pgate(di, cj);
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