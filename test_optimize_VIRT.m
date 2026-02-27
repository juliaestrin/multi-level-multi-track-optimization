% Test Script: optimize_VIRT.m
%
% Author:  Julia Estrin
% Date:    02-17-2026
% Updated: 02-27-2026
%
% Description:
%   Optimizes a 6.25 kW VIRT transformer design by sweeping over core loss
%   density (Pv_max) and window height. Finds the design that minimizes
%   total loss (core + copper) and displays the optimal parameters.
%   Includes resonant frequency calculation and 3D visualization.

clear; close all; clc;

fprintf('\n========================================\n');
fprintf('VIRT TRANSFORMER OPTIMIZATION\n');
fprintf('========================================\n');

%% Operating Point
Vo  = 48;                       % [V]   Output voltage per track
nt  = 2;                        % [-]   Number of secondary tracks
Lu  = 0.63e-6 * 5;              % [H]   Inductance requirement for 6.25 kW
Llk = 0.63e-6;                  % [H]   Leakage inductance
I   = 21.4 * sqrt(2);           % [A]   Peak primary current for 6.25 kW
fsw = 500e3;                    % [Hz]  Switching frequency
f   = 2 * fsw;                  % [Hz]  Transformer core frequency (2x due to FCML)
Np  = 4;                        % [-]   Number of primary turns (4:1/2 VIRT)

fprintf('\nOperating Point:\n');
fprintf('  Vo:             %.0f V\n', Vo);
fprintf('  nt:             %d tracks\n', nt);
fprintf('  Lu:             %.2f µH\n', Lu * 1e6);
fprintf('  Llk:            %.2f µH\n', Llk * 1e6);
fprintf('  I_peak:         %.2f A\n', I);
fprintf('  fsw:            %.0f kHz\n', fsw / 1e3);
fprintf('  f (core):       %.0f MHz\n', f / 1e6);
fprintf('  Np:             %d turns\n', Np);

%% Ferrite Material Selection
material_name = 'F80';        % Options: 'ML91S', 'F80', 'N87', 'N97', '3F4'
material = get_core_material(material_name);

k    = material.k;              % Steinmetz coefficient k   (P = k * B^beta)
beta = material.beta;           % Steinmetz exponent beta
uc   = material.uc;             % [-]  Relative permeability of core material

fprintf('\nCore Material:\n');
fprintf('  Material:       %s\n', material.name);
fprintf('  k:              %.4e\n', k);
fprintf('  beta:           %.4f\n', beta);
fprintf('  uc:             %d\n', uc);

%% Copper Parameters
rho_cu   = 2.2e-8;              % [Ohm·m] Resistivity at 100°C
sigma_cu = 1 / rho_cu;          % [S/m]   Conductivity at 100°C
u0       = 4 * pi * 1e-7;       % [H/m]   Permeability of free space

%% Winding Configuration
w_scale  = 1/2;                 % [-]     Winding width scale factor (width/length)
stackup  = '5layer';            % [-]     Winding stackup configuration

% Available stackup configurations:
%   '3layer'             - P-S-P
%   '5layer'             - P-P-P-P-S (non-interleaved)
%   '6layer'             - P-P-S-S-P-P (non-interleaved)
%   '6layer_interleaved' - P-S-P-P-S-P
%   '7layer_interleaved' - P-S-P-S-P-S-P
%   '8layer_interleaved' - P-S-P-S-P-S-P-S

fprintf('\nWinding Configuration:\n');
fprintf('  Stackup:        %s\n', stackup);
fprintf('  Width scale:    %.2f\n', w_scale);

%% Transformer Fixed Parameters
w_b      = 4e-3;                % [m]     Winding window breadth/depth
t_cu_pri = 2 * 35e-6;           % [m]     Primary copper thickness (2 oz)
t_cu_sec = 2 * 35e-6;           % [m]     Secondary copper thickness (2 oz)

fprintf('  w_b:            %.1f mm\n', w_b * 1e3);
fprintf('  t_cu_pri:       %.0f µm (%.1f oz)\n', t_cu_pri * 1e6, t_cu_pri / 35e-6);
fprintf('  t_cu_sec:       %.0f µm (%.1f oz)\n', t_cu_sec * 1e6, t_cu_sec / 35e-6);

%% Package Design Parameters
design_params = struct( ...
    'Vo',       Vo,       ...
    'nt',       nt,       ...
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
    'w_scale',  w_scale,  ...
    't_cu_pri', t_cu_pri, ...
    't_cu_sec', t_cu_sec, ...
    'stackup',  stackup   ...
);

%% Optimization Sweep
fprintf('\n========================================\n');
fprintf('RUNNING OPTIMIZATION SWEEP\n');
fprintf('========================================\n');

w_h_max       = 20e-3;                     % [m] Maximum window height
Pv_max_list   = linspace(50e3, 500e3, 50); % [W/m^3] Core loss density sweep
w_height_list = linspace(5e-3, w_h_max, 50); % [m]   Window height sweep

fprintf('  Pv_max range:   %.0f - %.0f kW/m³ (%d points)\n', ...
    min(Pv_max_list)/1e3, max(Pv_max_list)/1e3, length(Pv_max_list));
fprintf('  w_height range: %.1f - %.1f mm (%d points)\n', ...
    min(w_height_list)*1e3, max(w_height_list)*1e3, length(w_height_list));
fprintf('  Total designs:  %d\n', length(Pv_max_list) * length(w_height_list));

tic;
opt = optimize_VIRT(Pv_max_list, w_height_list, design_params);
elapsed = toc;

fprintf('  Optimization completed in %.2f seconds.\n', elapsed);

%% Parasitic Capacitance and Resonant Frequency
% Cps   = calculate_Cps_3Layer(opt.opt_design.l_winding, opt.opt_design.w_winding, opt.opt_design.w_core);
% f_res = 1 / (2 * pi * sqrt(Cps * Llk));

%% Display Optimal Design Results
fprintf('\n========================================\n');
fprintf('OPTIMAL DESIGN RESULTS\n');
fprintf('========================================\n');

fprintf('\n-- Operating Point --\n');
fprintf('  Pv_max_opt:     %.2f kW/m³\n', opt.Pv_max_opt / 1e3);
fprintf('  w_height_opt:   %.2f mm\n',    opt.w_height_opt * 1e3);
fprintf('  Bmax_opt:       %.4f T\n',     opt.Bmax_opt);

fprintf('\n-- Core Geometry --\n');
fprintf('  Ac:             %.2f mm²\n',   opt.opt_design.Ac * 1e6);
fprintf('  r_centerpost:   %.2f mm\n',    opt.opt_design.r_centerpost * 1e3);
fprintf('  w_core:         %.2f mm\n',    opt.opt_design.w_core * 1e3);
fprintf('  l_leg:          %.2f mm\n',    opt.opt_design.l_leg * 1e3);
fprintf('  lg (air gap):   %.4f mm\n',    opt.opt_design.lg * 1e3);
fprintf('  l_core:         %.2f mm\n',    opt.opt_design.l_core * 1e3);
fprintf('  h_core:         %.2f mm\n',    opt.opt_design.h_core * 1e3);
fprintf('  V_total:        %.2f cm³\n',   opt.opt_design.V_total * 1e6);
fprintf('  A_footprint:    %.2f cm²\n',   opt.opt_design.A_footprint * 1e4);

fprintf('\n-- Winding Geometry --\n');
fprintf('  l_winding:      %.2f mm\n',    opt.opt_design.l_winding * 1e3);
fprintf('  w_winding:      %.2f mm\n',    opt.opt_design.w_winding * 1e3);

fprintf('\n-- Resistance --\n');
fprintf('  Rdc_pri:        %.4f mΩ\n',    opt.opt_design.Rdc_pri * 1e3);
fprintf('  Rdc_sec:        %.4f mΩ\n',    opt.opt_design.Rdc_sec * 1e3);

fprintf('\n-- Losses --\n');
fprintf('  P_core:         %.4f W  (%.1f%%)\n', opt.P_core_min, 100*opt.P_core_min/opt.P_total_min);
fprintf('  P_copper:       %.4f W  (%.1f%%)\n', opt.P_copper_min, 100*opt.P_copper_min/opt.P_total_min);
fprintf('  P_total:        %.4f W\n', opt.P_total_min);

% fprintf('\n-- Parasitics --\n');
% fprintf('  Cps:            %.4f nF\n',    Cps * 1e9);
% fprintf('  f_res (Llk-Cps):%.4f MHz\n',   f_res * 1e-6);

fprintf('\n========================================\n\n');

%% 3D Visualization of Optimal Design
fprintf('Generating 3D visualization...\n');
core3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt);
fprintf('3D visualization complete.\n');
fprintf('========================================\n\n');