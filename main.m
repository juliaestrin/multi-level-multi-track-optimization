%% ============= Add Path =============
addpath('Switch Functions');
addpath('Transformer Design');
addpath('LLC Design');


%% Design Specifications
Vin_nom     = 1500;     % [V]    Nominal input voltage
Vo_nom      = 48;       % [V]    Nominal output voltage
nt          = 2;        % [-]    Number of tracks

% Select FCML switching freqency
fsw         = 500e3;    % [Hz]   FCML switching frequency
f0 = 2 * fsw; %[Hz] LLC and secondary side frequency 

% Select Converter Power Level 
Pmax        = 6.25e3;   % [W]    Maximum output power
Pmin        = 0.1*Pmax; % [W]    Minimum output power (10% of rated)

% Transformer Design Specifications 
material_name = 'F80';        % Options: 'ML91S', 'F80', '3F46'
w_h_max = 100e-3;               % [m] max window height (width) 
w_scale  = 1/2;                 % [-]     Winding width scale factor (width/length)
stackup = '5layer'; 
%   Supported stackup configurations:
%     '3layer'             - 3-layer: P-S-P
%     '5layer'             - 5-layer: P-P-P-P-S
%     '5layer_interleaved' - 5-layer: P-P-S-P-P
%     '6layer'             - 6-layer non-interleaved: P-P-S-S-P-P
%     '6layer_interleaved' - 6-layer interleaved: P-S-P-P-S-P
%     '7layer_interleaved' - 7-layer interleaved: P-S-P-S-P-S-P
%     '8layer_interleaved' - 8-layer interleaved: P-S-P-S-P-S-P-S

% Transformer Fixed Parameters
w_b      = 4e-3;                % [m]     Winding window breadth/depth
t_cu_pri = 2 * 35e-6;           % [m]     Primary copper thickness (2 oz)
t_cu_sec = 2 * 35e-6;           % [m]     Secondary copper thickness (2 oz)

% LLC Design Specifications 
Mg_nom      = 1.0;      % [-]    Nominal LLC gain (unity at resonance)
percentReg  = 0.1;      % [-]    Line regulation ±10%
f_per       = 0.25;     % [-]    Frequency range ±25%
Ln          = 5;        % [-]    Inductance ratio Lm/Lr

% Constants 
rho_cu   = 2.2e-8;              % [Ohm·m] Resistivity at 100°C
sigma_cu = 1 / rho_cu;          % [S/m]   Conductivity at 100°C
u0       = 4 * pi * 1e-7;       % [H/m]   Permeability of free space
%% Design LLC 
LLC_design = designLLC(Vin_nom, Vo_nom, Mg_nom, nt, percentReg, fsw, f_per, Pmax, Pmin, Ln);

% Extract results 
Lu = LLC_design.Lm;  % [H]  Inductance requirement for 6.25 kW
Llk = LLC_design.Lr; 
Ir_rms = LLC_design.Ir_rms; 
Ir_pk = Ir_rms * sqrt(2);
N = LLC_design.N; % total primary side turns for nt tracks; 
Np = N / nt; 


%% Design Transformer 

material = get_core_material(material_name);

k    = material.k;              % Steinmetz coefficient k   (P = k * B^beta)
beta = material.beta;           % Steinmetz exponent beta
uc   = material.uc;             % [-]  Relative permeability of core material

% Package Design Parameters
design_params = struct( ...
    'Vo',       Vo_nom,       ...
    'nt',       nt,       ...
    'Lu',       Lu,       ...
    'I',        Ir_pk,        ...
    'f',        f0,        ...
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

% Optimization Sweep
Pv_max_list   = linspace(50e3, 500e3, 500);     % [W/m^3] core loss density sweep
w_height_list = linspace(5e-3, w_h_max, 50);     % [m]     window height sweep
opt = optimize_VIRT(Pv_max_list, w_height_list, design_params);

core3Dfigure(opt.opt_design, opt.Pv_max_opt, opt.w_height_opt);

% Cps = calculcate_Cps_3Layer(opt.opt_design.l_winding, opt.opt_design.w_winding, opt.opt_design.w_core); 
% f_res = 1/(2*pi*sqrt(Cps * Llk));

%% ============= Call Functions =============
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
out1 = analyzePriSwitches(Pmax, fsw, Ir_rms, 2, 10, ...
    [1 2 3 4 5 6 7 8], [1 2 3 4 5 6 7 8], 500, ...
    'GaN Data.xlsx', 'SiC Data.xlsx');

% Secondary Side
out2 = analyzeSecSwitches(Pmax, f0, 1, 10, [4 6 8], [4 6]);

%%
% Efficiency Calculation
% eff = calcEfficiency(out1, out2, "pareto", "minLoss", 6.25e3, opt.P_total_min, opt.opt_design.A_footprint*1e6);