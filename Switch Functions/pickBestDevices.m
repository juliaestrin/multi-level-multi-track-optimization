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