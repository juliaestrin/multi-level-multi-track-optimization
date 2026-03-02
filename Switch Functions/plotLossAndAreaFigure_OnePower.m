function plotLossAndAreaFigure_OnePower(C, rank, marker_by_jj, figTitle, figNum)
    figure(figNum); clf;
    tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

    jj_set = C.jj_set;
    jj_list = rank.jj_list;
    jj_list = jj_list(ismember(jj_list, jj_set));

    nexttile; hold on; grid on;
    for i = 1:C.N_total

        fillFlag = true;
        if isfield(C,'MarkerFilled') && numel(C.MarkerFilled) >= i
            fillFlag = C.MarkerFilled(i);
        end

        if fillFlag
            ls = '-';
            mk = 'o';
            mfc = C.colors(i,:);
        else
            ls = ':';
            mk = 's';
            mfc = 'none';
        end

        plot(jj_set, C.P_loss(i,:), ...
            'LineStyle', ls, 'Marker', mk, ...
            'Color', C.colors(i,:), ...
            'MarkerFaceColor', mfc, ...
            'MarkerEdgeColor', C.colors(i,:), ...
            'LineWidth', 1.5, 'MarkerSize', 6);
    end
    xlabel('Number of Parallel Switches','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);
    legend(C.PartNumber,'Location','best','FontSize',10);

    nexttile; hold on; grid on;
    msize = 90;
    h_jj = gobjects(numel(jj_list),1);

    for k = 1:numel(jj_list)
        jj_abs = jj_list(k);
        col = find(jj_set == jj_abs, 1);

        mk = 'o';
        if jj_abs <= numel(marker_by_jj) && ~isempty(marker_by_jj{jj_abs})
            mk = marker_by_jj{jj_abs};
        end

        for i = 1:C.N_total
            if isnan(C.P_loss(i,col)), continue; end

            fillFlag = true;
            if isfield(C,'MarkerFilled') && numel(C.MarkerFilled) >= i
                fillFlag = C.MarkerFilled(i);
            end

            if fillFlag
                hh = scatter(C.Area(i,col), C.P_loss(i,col), msize, ...
                    'Marker', mk, ...
                    'MarkerFaceColor', C.colors(i,:), ...
                    'MarkerEdgeColor', C.colors(i,:), ...
                    'LineWidth', 1.2);
            else
                hh = scatter(C.Area(i,col), C.P_loss(i,col), msize, ...
                    'Marker', mk, ...
                    'MarkerFaceColor', 'none', ...
                    'MarkerEdgeColor', C.colors(i,:), ...
                    'LineWidth', 1.2);
            end

            if ~isgraphics(h_jj(k)), h_jj(k) = hh; end
        end
    end

    xlabel('Total Area [mm$^2$]','Interpreter','latex','FontSize',15);
    ylabel('Total Power Loss [W]','Interpreter','latex','FontSize',15);

    if ~isempty(jj_list)
        ok = isgraphics(h_jj);
        legend(h_jj(ok), compose('%d-parallel', jj_list(ok)), ...
            'Location','best','Orientation','horizontal','FontSize',10);
    end

    sgtitle(figTitle,'Interpreter','latex','FontSize',15);
end