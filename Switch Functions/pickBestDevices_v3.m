% Comparing to v2:
% Supports both:
%   1) scalar jj candidates, e.g. jj_set = [1 2 3 4]
%   2) pair candidates, e.g. jj_set = [jj25 jj75] as an N x 2 matrix
%
% Inputs:
%   AreaFull, LossFull, LossAreaFull, PcondFull, PoffFull, PgateFull
%       : size = [nDevices x nCandidates]
%
%   pri_sw_count
%       : for scalar-jj topologies, per-position loss is multiplied by pri_sw_count
%         to get total primary-side loss
%       : for pair-jj topologies, pass [] or 1 if PcondFull/PoffFull/PgateFull
%         are already total primary-side values
%
%   jj_set
%       : either 1 x N / N x 1 vector   OR   N x 2 matrix
%
%   compare_list
%       : subset of jj_set to compare
%         - scalar case: vector, e.g. [2 4 6]
%         - pair case:   M x 2 matrix, e.g. [1 2; 2 2; 2 4]
%
function best = pickBestDevices_v3(AreaFull, LossFull, LossAreaFull, ...
    PcondFull, PoffFull, PgateFull, pri_sw_count, jj_set, compare_list, Meta)

    %% ===================== Detect candidate format =====================
    isPairMode = ismatrix(jj_set) && size(jj_set,2) == 2;

    if isPairMode
        if isempty(compare_list)
            compare_list = jj_set;
        end
        if size(compare_list,2) ~= 2
            error('For pair mode, compare_list must be an N x 2 matrix of [jj25, jj75].');
        end
    else
        if isempty(compare_list)
            compare_list = jj_set;
        end
        compare_list = compare_list(:).';   % force row vector
    end

    %% ===================== Find selected columns =====================
    if isPairMode
        % jj_set is N x 2, compare_list is M x 2
        cols = zeros(1, size(compare_list,1));
        for k = 1:size(compare_list,1)
            tf = ismember(jj_set, compare_list(k,:), 'rows');
            idx = find(tf, 1);
            if isempty(idx)
                error('compare_list row [%g %g] was not found in jj_set.', ...
                    compare_list(k,1), compare_list(k,2));
            end
            cols(k) = idx;
        end
    else
        cols = zeros(1, numel(compare_list));
        for k = 1:numel(compare_list)
            idx = find(jj_set == compare_list(k), 1);
            if isempty(idx)
                error('compare_list entry %g was not found in jj_set.', compare_list(k));
            end
            cols(k) = idx;
        end
    end

    %% ===================== Slice selected columns =====================
    Loss     = LossFull(:, cols);
    Area     = AreaFull(:, cols);
    LossArea = LossAreaFull(:, cols);

    % Keep breakdown consistent with Loss
    if isPairMode
        % In pair mode, PcondFull/PoffFull/PgateFull should already be TOTAL
        % primary-side values for that [jj25, jj75] design point.
        Pcond = PcondFull(:, cols);
        Poff  = PoffFull(:, cols);
        Pgate = PgateFull(:, cols);
    else
        % Original scalar-jj behavior
        Pcond = pri_sw_count * PcondFull(:, cols);
        Poff  = pri_sw_count * PoffFull(:, cols);
        Pgate = pri_sw_count * PgateFull(:, cols);
    end

    %% ===================== Feasibility mask =====================
    feas = ~isnan(Loss);

    Pcond(~feas) = NaN;
    Poff(~feas)  = NaN;
    Pgate(~feas) = NaN;

    %% ===================== Minimum loss =====================
    Loss1 = Loss;
    Loss1(~feas) = NaN;
    [~, lin] = min(Loss1(:), [], 'omitnan');

    if isempty(lin) || all(isnan(Loss1(:)))
        best.minLoss = struct([]);
        best.minArea = struct([]);
        best.minLossArea = struct([]);
        best.pareto = struct([]);
        best.compare_list = compare_list;
        return;
    end

    [dev, colLocal] = ind2sub(size(Loss1), lin);
    best.minLoss = packOneLocal_v3(dev, colLocal, compare_list, ...
        Loss, Area, Pcond, Poff, Pgate, Meta, isPairMode);

    %% ===================== Minimum area =====================
    Area1 = Area;
    Area1(~feas) = NaN;
    [~, lin] = min(Area1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(Area1), lin);
    best.minArea = packOneLocal_v3(dev, colLocal, compare_list, ...
        Loss, Area, Pcond, Poff, Pgate, Meta, isPairMode);

    %% ===================== Minimum loss-area product =====================
    LA1 = LossArea;
    LA1(~feas) = NaN;
    [~, lin] = min(LA1(:), [], 'omitnan');
    [dev, colLocal] = ind2sub(size(LA1), lin);
    best.minLossArea = packOneLocal_v3(dev, colLocal, compare_list, ...
        Loss, Area, Pcond, Poff, Pgate, Meta, isPairMode);

    %% ===================== Pareto front =====================
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

    if isPairMode
        paretoList = repmat(struct( ...
            'device_index', [], ...
            'device_name', "", ...
            'jj25', [], ...
            'jj75', [], ...
            'loss_W', [], ...
            'area_mm2', [], ...
            'loss_area_W_mm2_product', [], ...
            'P_cond_W', [], ...
            'P_off_W', [], ...
            'P_gate_W', []), ...
            numel(pDevIdx), 1);
    else
        paretoList = repmat(struct( ...
            'device_index', [], ...
            'device_name', "", ...
            'jj', [], ...
            'loss_W', [], ...
            'area_mm2', [], ...
            'loss_area_W_mm2_product', [], ...
            'P_cond_W', [], ...
            'P_off_W', [], ...
            'P_gate_W', []), ...
            numel(pDevIdx), 1);
    end

    for m = 1:numel(pDevIdx)
        di = pDevIdx(m);
        cj = pColIdx(m);

        paretoList(m).device_index = di;
        paretoList(m).device_name  = string(Meta.Name(di));

        if isPairMode
            paretoList(m).jj25 = compare_list(cj,1);
            paretoList(m).jj75 = compare_list(cj,2);
        else
            paretoList(m).jj = compare_list(cj);
        end

        paretoList(m).loss_W  = Loss(di, cj);
        paretoList(m).area_mm2 = Area(di, cj);
        paretoList(m).loss_area_W_mm2_product = Area(di, cj) * Loss(di, cj);

        paretoList(m).P_cond_W = Pcond(di, cj);
        paretoList(m).P_off_W  = Poff(di, cj);
        paretoList(m).P_gate_W = Pgate(di, cj);
    end

    %% ===================== Sort Pareto list =====================
    if ~isempty(paretoList)
        A = [paretoList.area_mm2].';
        L = [paretoList.loss_W].';

        if isPairMode
            J25 = [paretoList.jj25].';
            J75 = [paretoList.jj75].';
            [~, ord] = sortrows([A L J25 J75], [1 2 3 4]);
        else
            J = [paretoList.jj].';
            [~, ord] = sortrows([A L J], [1 2 3]);
        end

        paretoList = paretoList(ord);
    end

    %% ===================== Outputs =====================
    best.pareto = paretoList;
    best.compare_list = compare_list;
    best.isPairMode = isPairMode;
end