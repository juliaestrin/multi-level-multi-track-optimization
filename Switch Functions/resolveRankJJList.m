function jj_list_out = resolveRankJJList(mode, jj_list_in, max_para, selected_para)
    jj_list_in = unique(jj_list_in(:).','stable');

    if mode == 1
        jj_list_out = jj_list_in(jj_list_in >= 1 & jj_list_in <= max_para);
        if isempty(jj_list_out), jj_list_out = 1:max_para; end
    elseif mode == 2
        sel = unique(selected_para(:).','stable');
        overlap = jj_list_in(ismember(jj_list_in, sel));
        if isempty(overlap), jj_list_out = sel; else, jj_list_out = overlap; end
    else
        error('mode must be 1 or 2');
    end
end