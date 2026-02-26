function isPareto = computeParetoFront(area_vec, loss_vec)
    n = numel(area_vec);
    isPareto = true(n,1);
    for i = 1:n
        if ~isPareto(i), continue; end
        ai = area_vec(i);
        li = loss_vec(i);
        dominates_i = (area_vec <= ai) & (loss_vec <= li) & ((area_vec < ai) | (loss_vec < li));
        if any(dominates_i)
            isPareto(i) = false;
        end
    end
end