out1 = analyzePriSwitches_v3("Multilevel Multitrack", 6.25e3, 500e3, 21.4142, 1, 8, ...
    [], [], 10000, ...
    'GaN Data tf.xlsx', 'SiC Data tf.xlsx');
%
% Secondary Side
out2 = analyzeSecSwitches(6.25e3, 1000e3, 1, 10, [], [4 6 8]);

eff = calcEfficiency_v4(out1, out2, "pareto", "pareto", 6.25e3, 0, 0);

