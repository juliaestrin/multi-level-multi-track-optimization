% out1 = analyzePriSwitches_v3("Multilevel Multitrack", 6.25e3, 500e3, 21.4142, 1, 8, ...
%     [], [], 10000, ...
%     'GaN Data tf.xlsx', 'SiC Data tf.xlsx');
out1 = analyzePriSwitches_v3("Multilevel Multitrack", 6.25e3, 500e3, 21.447071, 1, 8, ...
    [], [], 10000, ...
    'GaN Data tf.xlsx', 'SiC Data tf.xlsx');

try
    pareto = out1.best.pareto;
    T_pareto = struct2table(pareto);

    if topology == "Multilevel Multitrack"
        wantedVars = {'device_name','jj25','jj75','P_cond_W','P_off_W','P_gate_W','loss_W','area_mm2'};
        hasVars = ismember(wantedVars, T_pareto.Properties.VariableNames);

        if all(hasVars)
            T_pareto = T_pareto(:, wantedVars);
            T_pareto.Properties.VariableNames = { ...
                'DeviceName', 'jj25', 'jj75', 'P_cond', 'P_off', ...
                'P_gate', 'TotalLoss', 'TotalArea_mm2'};
        end
    else
        wantedVars = {'device_name','jj','P_cond_W','P_off_W','P_gate_W','loss_W','area_mm2'};
        hasVars = ismember(wantedVars, T_pareto.Properties.VariableNames);

        if all(hasVars)
            T_pareto = T_pareto(:, wantedVars);
            T_pareto.Properties.VariableNames = { ...
                'DeviceName', 'ParallelCount', 'P_cond', 'P_off', ...
                'P_gate', 'TotalLoss', 'TotalArea_mm2'};
        end
    end

    disp(T_pareto);
catch
    warning('Could not display primary pareto table for %s.', topology);
end

% Secondary Side
out2 = analyzeSecSwitches(6.25e3, 1000e3, 1, 10, [], [4 6 8]);

eff = calcEfficiency_v4(out1, out2, "pareto", "pareto", 6.25e3, 0, 0);

