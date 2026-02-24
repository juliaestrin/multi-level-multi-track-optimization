% Test Script: optimize_VIRT.m
%
% Author:  Julia Estrin
% Date:    02-17-2026
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

%% Optimization Sweep
Pv_max_list   = linspace(50e3, 500e3, 50);     % [W/m^3] core loss density sweep
w_height_list = linspace(5e-3, 5*w_b, 50);     % [m]     window height sweep
opt = optimize_VIRT(Pv_max_list, w_height_list, design_params);

Cps = calculcate_Cps_3Layer(opt.opt_design.l_winding, opt.opt_design.w_winding, opt.opt_design.w_core); 
f_res = 1/(2*pi*sqrt(Cps * Lu));

%% Display Optimal Design Results
fprintf('\n===== Optimal Design Results =====\n');
fprintf('  Pv_max_opt:     %.2f kW/m^3\n', opt.Pv_max_opt   / 1e3);
fprintf('  w_height_opt:   %.2f mm\n',      opt.w_height_opt * 1e3);
fprintf('  Bmax_opt:       %.4f T\n',        opt.Bmax_opt);
fprintf('  P_total_min:    %.4f W\n',        opt.P_total_min);
fprintf('  P_core_min:     %.4f W\n',        opt.P_core_min);
fprintf('  P_copper_min:   %.4f W\n',        opt.P_copper_min);
fprintf('  f_res:          %.4f MHz\n',        f_res*1e-6);
fprintf('==================================\n\n');


%% 3D Visualization of Optimal Design
core3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt);