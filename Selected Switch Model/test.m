clear; clc; close all;

data = readtable('plot-data.csv');

x = data{:,1};   % Tj [degC]
y = data{:,2};   % normalized Rds(on)

% Sort data by x
% [x, idx] = sort(x);
% y = y(idx);

% 4th-order polynomial fitting
p = polyfit(x, y, 4);

% Generate smooth fitted curve
xfit = linspace(min(x), max(x), 500);
yfit = polyval(p, xfit);

% Plot
figure;
plot(x, y, 'o', 'LineWidth', 1.5);
hold on;
plot(xfit, yfit, 'r-', 'LineWidth', 2);
grid on;

xlabel('$T_j$ [$^\circ$C]', 'Interpreter', 'latex');
ylabel('$R_{DS(on)}$ [normalized]', 'Interpreter', 'latex');
legend('Extracted data', '4th-order fit', 'Location', 'northwest');

% Display equation coefficients
fprintf('y = %.4e*x^4 + %.4e*x^3 + %.4e*x^2 + %.4e*x + %.4e\n', ...
    p(1), p(2), p(3), p(4), p(5));

% Calculate R-square
y_pred = polyval(p, x);
SS_res = sum((y - y_pred).^2);
SS_tot = sum((y - mean(y)).^2);
R2 = 1 - SS_res/SS_tot;

fprintf('R^2 = %.6f\n', R2);