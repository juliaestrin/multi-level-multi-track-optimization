function s = packOneLocal(devIdx, jj_abs, compare_list, LossSub, AreaSub, Meta)
    colLocal = find(compare_list == jj_abs, 1);
    if isempty(colLocal)
        error('packOneLocal: jj_abs=%d not found in compare_list=[%s].', jj_abs, num2str(compare_list));
    end

    s = struct();
    s.device_index = devIdx;
    s.device_name  = string(Meta.Name(devIdx));
    s.jj           = jj_abs;

    s.loss_W  = LossSub(devIdx, colLocal);
    s.area_mm2 = AreaSub(devIdx, colLocal);
    s.loss_area_W_mm2_product = s.loss_W * s.area_mm2;
end