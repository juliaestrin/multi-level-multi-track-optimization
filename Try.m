addpath('Switch Functions');

% Primary Side: 
% 1) GaN only:
% out1 = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
%     [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
%     'GaN Data.xlsx', []);
%
% 2) SiC only:
% out1 = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
%     [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
%     [], 'SiC Data.xlsx');
%
% 3) Both:
out1 = analyzePriSwitches(6.25e3, 1000e3, 21.4142, 2, 10, ...
    [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
    'GaN Data.xlsx', 'SiC Data.xlsx');
%
% Secondary Side
out2 = analyzeSecSwitches(6.25e3, 1000e3, 1, 10, [4 6 8], [4 6]);

%% Operating Point
Lu  = 0.63e-6 * 5;             % [H]  Inductance requirement for 6.25 kW
Llk = 0.63e-6; 
I   = 21.4 * sqrt(2);          % [A]  Peak primary current for 6.25 kW
fsw = 500e3;                   % [Hz] Switching frequency
f   = 2 * fsw;                 % [Hz] Transformer core frequency (2x switching due to FCML frequency doubling)
Np  = 4;                       % [-]  Number of primary turns (4 : 1/2 VIRT)
w_h_max = 100e-3;               % [m] max window height (width) 

%% Ferrite Material Parameters (Proterials ML91S)
k    = 1606769230.46174;        % Steinmetz coefficient k   (P = k * B^beta)
beta = 3.19564613808786;        % Steinmetz exponent beta
uc   = 900;                     % [-]  Relative permeability of core material

%% Copper Parameters
rho_cu   = 2.2e-8;              % [Ohm·m] Resistivity at 100°C
sigma_cu = 1 / rho_cu;          % [S/m]   Conductivity at 100°C
u0       = 4 * pi * 1e-7;       % [H/m]   Permeability of free space

stackup = '5layer'; 
%   Supported stackup configurations:
%     '3layer'             - 3-layer: P-S-P
%     '5layer'             - 5-layer: P-P-P-P-S
%     '5layer_interleaved' - 5-layer: P-P-S-P-P
%     '6layer'             - 6-layer non-interleaved: P-P-S-S-P-P
%     '6layer_interleaved' - 6-layer interleaved: P-S-P-P-S-P
%     '7layer_interleaved' - 7-layer interleaved: P-S-P-S-P-S-P
%     '8layer_interleaved' - 8-layer interleaved: P-S-P-S-P-S-P-S

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
    't_cu_sec', t_cu_sec,  ...
    'stackup', stackup ...
);

%% Optimization Sweep
Pv_max_list   = linspace(50e3, 500e3, 500);     % [W/m^3] core loss density sweep
w_height_list = linspace(5e-3, w_h_max, 50);     % [m]     window height sweep
opt = optimize_VIRT(Pv_max_list, w_height_list, design_params);

Cps = calculcate_Cps_3Layer(opt.opt_design.l_winding, opt.opt_design.w_winding, opt.opt_design.w_core); 
f_res = 1/(2*pi*sqrt(Cps * Llk));

% Efficiency Calculation
eff = calcEfficiency(out1, out2, "pareto", "minLoss", 6.25e3, opt.P_total_min, opt.opt_design.A_footprint*(10^6));