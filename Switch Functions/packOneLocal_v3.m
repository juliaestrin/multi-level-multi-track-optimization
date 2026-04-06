function s = packOneLocal_v3(devIdx, colLocal, compare_list, LossSub, AreaSub, PcondSub, PoffSub, PgateSub, Meta, isPairMode)

    % Basic index checks
    if colLocal < 1 || colLocal > size(LossSub,2)
        error('packOneLocal_v3: colLocal=%d is out of range.', colLocal);
    end

    if devIdx < 1 || devIdx > size(LossSub,1)
        error('packOneLocal_v3: devIdx=%d is out of range.', devIdx);
    end

    s = struct();
    s.device_index = devIdx;
    s.device_name  = string(Meta.Name(devIdx));

    if isPairMode
        if size(compare_list,2) ~= 2
            error('packOneLocal_v3: compare_list must be N x 2 in pair mode.');
        end

        s.jj25 = compare_list(colLocal, 1);
        s.jj75 = compare_list(colLocal, 2);
    else
        compare_list = compare_list(:).';
        s.jj = compare_list(colLocal);
    end

    s.loss_W   = LossSub(devIdx, colLocal);
    s.area_mm2 = AreaSub(devIdx, colLocal);
    s.loss_area_W_mm2_product = s.loss_W * s.area_mm2;

    s.P_cond_W = PcondSub(devIdx, colLocal);
    s.P_off_W  = PoffSub(devIdx, colLocal);
    s.P_gate_W = PgateSub(devIdx, colLocal);
end