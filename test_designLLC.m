% Test Script: designLLC.m
%
% Author:  Julia Estrin
% Date:    02-27-2026
%
% Description:
%   Tests the LLC resonant converter design function for an FCML VIRT
%   topology. Designs the resonant tank for a 6.25 kW converter with
%   1500V nominal input and 48V output.

clear; close all; clc;

fprintf('\n========================================\n');
fprintf('LLC RESONANT CONVERTER DESIGN\n');
fprintf('========================================\n');

%% Design Specifications
Mg_nom      = 1.0;      % [-]    Nominal LLC gain (unity at resonance)
nt          = 2;        % [-]    Number of secondary tracks
Vin_nom     = 1500;     % [V]    Nominal input voltage
Vo_nom      = 48;       % [V]    Nominal output voltage
percentReg  = 0.1;      % [-]    Line regulation ±10%
Pmax        = 6.25e3;   % [W]    Maximum output power
Pmin        = 0.1*Pmax; % [W]    Minimum output power (10% of rated)
fsw         = 500e3;    % [Hz]   FCML switching frequency
Ln          = 5;        % [-]    Inductance ratio Lm/Lr
f_per       = 0.25;     % [-]    Frequency range ±25%

fprintf('\nInput Specifications:\n');
fprintf('  Vin_nom:      %.0f V\n', Vin_nom);
fprintf('  Vo_nom:       %.0f V\n', Vo_nom);
fprintf('  Pmax:         %.2f kW\n', Pmax / 1e3);
fprintf('  Pmin:         %.2f kW\n', Pmin / 1e3);
fprintf('  fsw:          %.0f kHz\n', fsw / 1e3);
fprintf('  Ln (Lm/Lr):   %.1f\n', Ln);
fprintf('  Reg:          ±%.0f%%\n', percentReg * 100);
fprintf('  Freq range:   ±%.0f%%\n', f_per * 100);

%% Run LLC Design
result = designLLC(Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, f_per, Pmax, Pmin, Ln);

%% Display Results
fprintf('\n========================================\n');
fprintf('DESIGN RESULTS\n');
fprintf('========================================\n');

fprintf('\n-- Transformer --\n');
fprintf('  Turns ratio N:   %d:1 (total)\n', result.N);
fprintf('  Tracks nt:       %d\n', nt);

fprintf('\n-- Operating Frequencies --\n');
fprintf('  f0 (resonant):   %.2f MHz\n', result.f0 / 1e6);
fprintf('  fn_min:          %.3f  (%.2f MHz)\n', result.fn_min, result.fn_min * result.f0 / 1e6);
fprintf('  fn_max:          %.3f  (%.2f MHz)\n', result.fn_max, result.fn_max * result.f0 / 1e6);

fprintf('\n-- Gain Range --\n');
fprintf('  Mg_min:          %.4f\n', result.Mg_min);
fprintf('  Mg_max:          %.4f\n', result.Mg_max);

fprintf('\n-- Resonant Tank Components --\n');
fprintf('  Lr:              %.4f µH\n', result.Lr * 1e6);
fprintf('  Cr:              %.4f µF\n', result.Cr * 1e6);
fprintf('  Lm:              %.4f µH\n', result.Lm * 1e6);
fprintf('  Qe_max:          %.4f\n', result.Qe_max);

fprintf('\n-- Equivalent Resistance --\n');
fprintf('  Re_min:          %.4f Ohm  (at Pmax)\n', result.Re_min);
fprintf('  Re_max:          %.4f Ohm  (at Pmin)\n', result.Re_max);

fprintf('\n-- RMS Current --\n');
fprintf('  Ir_rms:          %.2f A  (primary, worst case)\n', result.Ir_rms);

fprintf('\n========================================\n\n');
