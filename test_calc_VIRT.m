% Test Script: calculate_VIRT_design.m
%
% Author:  Julia Estrin
% Date:    02-12-2026
%
% Description:
%   Defines design parameters for a 6.25 kW VIRT transformer and calls
%   calculate_VIRT_design to evaluate core geometry, core loss, and AC
%   copper loss. Results are printed to the command window.

%% Operating Point
Lu  = 0.63e-6 * 5;             % [H]  Inductance requirement for 6.25 kW
I   = 21.4 * sqrt(2);          % [A]  Peak primary current for 6.25 kW
fsw = 500e3;                    % [Hz] Switching frequency
f   = 2 * fsw;                  % [Hz] Transformer core frequency (2x switching due to FCML frequency doubling)
Np  = 4;                        % [-]  Number of primary turns (4 : 1/2 VIRT)

%% Ferrite Material Parameters (Proterials ML91S)
k    = 1606769230.46174;        % Steinmetz coefficient k   (P = k * B^beta)
beta = 3.19564613808786;        % Steinmetz exponent beta
uc   = 900;                     % [-]  Relative permeability of core material

%% Copper Parameters
rho_cu   = 2.2e-8;              % [Ohm·m] Resistivity at 100°C
sigma_cu = 1 / rho_cu;          % [S/m]   Conductivity at 100°C
u0       = 4 * pi * 1e-7;       % [H/m]   Permeability of free space

%% Transformer Fixed Parameters
Pv_max   = 100e3;               % [W/m^3] Maximum allowable volumetric core loss
w_height = 25e-3;               % [m]     Winding window height
w_b      = 4e-3;                % [m]     Winding window breadth/depth
t_cu_pri = 2 * 35e-6;           % [m]     Primary copper thickness (2 oz)
t_cu_sec = 2 * 35e-6;           % [m]     Secondary copper thickness (2 oz)

%% Package Design Parameters
design_params = struct( ...
    'Lu',       Lu,       ...
    'I',        I,        ...
    'f',        f,        ...
    'Np',       Np,       ...
    'k',        k,        ...
    'beta',     beta,     ...
    'uc',       uc,       ...
    'rho_cu',   rho_cu,   ...
    'sigma_cu', sigma_cu, ...
    'u0',       u0,       ...
    'w_b',      w_b,      ...
    't_cu_pri', t_cu_pri, ...
    't_cu_sec', t_cu_sec  ...
);

%% Run Design
result = calculate_VIRT_design(Pv_max, w_height, design_params);

%% Display Results
fprintf('\n===== VIRT Design Results =====\n');
fprintf('\n-- Flux Density --\n');
fprintf('  Bmax:           %.4f T\n',   result.Bmax);

fprintf('\n-- Core Geometry --\n');
fprintf('  Ac:             %.4f mm^2\n', result.Ac           * 1e6);
fprintf('  r_centerpost:   %.4f mm\n',  result.r_centerpost * 1e3);
fprintf('  w_core:         %.4f mm\n',  result.w_core       * 1e3);
fprintf('  l_leg:          %.4f mm\n',  result.l_leg        * 1e3);
fprintf('  l_core:         %.4f mm\n',  result.l_core       * 1e3);
fprintf('  h_core:         %.4f mm\n',  result.h_core       * 1e3);
fprintf('  V_total:        %.4f cm^3\n', result.V_total     * 1e6);
fprintf('  A_footprint:    %.4f cm^2\n', result.A_footprint * 1e4);

fprintf('\n-- Winding Geometry --\n');
fprintf('  l_winding:      %.4f mm\n',  result.l_winding * 1e3);
fprintf('  w_winding:      %.4f mm\n',  result.w_winding * 1e3);

fprintf('\n-- Resistance --\n');
fprintf('  Rdc_pri:        %.4f mOhm\n', result.Rdc_pri * 1e3);
fprintf('  Rdc_sec:        %.4f mOhm\n', result.Rdc_sec * 1e3);

fprintf('\n-- Losses --\n');
fprintf('  P_core:         %.4f W\n',   result.P_core);
fprintf('  P_pri:          %.4f W\n',   result.P_pri);
fprintf('  P_sec:          %.4f W\n',   result.P_sec);
fprintf('  P_total:        %.4f W\n',   result.P_total);
fprintf('================================\n\n');

%% 3D Visualization
core3Dfigure(result, Pv_max, w_height);