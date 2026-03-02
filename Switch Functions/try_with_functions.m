% Primary Side
out1 = analyzePriSwitches(6.25e3, 500e3, 21.4142, 2, 10, ...
    [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
    'GaN Data.xlsx', 'SiC Data.xlsx');
%
% Secondary Side
out2 = analyzeSecSwitches(6.25e3, 1000e3, 1, 10, [4 6 8], [4 6]);


eff = calcEfficiency(out1, out2, "pareto", "minLoss", 6.25e3, 20, 0);