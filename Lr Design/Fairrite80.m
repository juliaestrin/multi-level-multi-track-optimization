B_data  = [0.0163351, 0.0204166, 0.023316, 0.0273354, 0.0318982, ...
    0.0372068, 0.0440032, 0.05062, 0.0594484];          % Tesla

Pv_data = [10.032, 20.931, 32.427, 54.861, 91.476, ...
    151.956, 263.914, 418.860, 715.316] * 1e3;          % W/m³

% Log-log fit: Pv = k * B^beta
% log(Pv) = log(k) + beta*log(B)
coeffs = polyfit(log(B_data), log(Pv_data), 1);

beta = coeffs(1)
k    = exp(coeffs(2))

% Verify
Pv_fit = k * B_data.^beta;
figure;
loglog(B_data, Pv_data, 'ro', B_data, Pv_fit, 'b-');
xlabel('B (T)'); ylabel('Pv (W/m^3)');
legend('Data', 'Fit'); grid on;
fprintf('k    = %.4e\n', k);
fprintf('beta = %.4f\n', beta);