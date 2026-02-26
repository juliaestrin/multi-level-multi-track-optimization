function C = buildSelectedTopN(Area, P_total_plot, jj_set, Meta, colors, rank)
    col_rank = find(jj_set == rank.jj_rank, 1);
    if isempty(col_rank)
        error('rank.jj_rank=%d is not in jj_set=[%s].', rank.jj_rank, num2str(jj_set));
    end

    loss_ref = P_total_plot(:, col_rank);
    valid = ~isnan(loss_ref);

    [~, idxs] = sort(loss_ref(valid), 'ascend');
    all_idx = find(valid);

    N_keep = min(rank.N_keep, numel(idxs));
    sel = all_idx(idxs(1:N_keep));

    C.Area       = Area(sel,:);
    C.P_loss     = P_total_plot(sel,:);
    C.PartNumber = string(Meta.Name(sel));
    C.N_total    = N_keep;
    C.jj_set     = jj_set;

    if size(colors,1) < numel(sel)
        C.colors = lines(numel(sel));
    else
        C.colors = colors(sel,:);
    end

    % ---- marker fill by tech if available ----
    if ismember('Tech', Meta.Properties.VariableNames)
        tech_sel = string(Meta.Tech(sel));
        C.MarkerFilled = (tech_sel == "GaN");   % GaN filled, SiC empty
    end
end
